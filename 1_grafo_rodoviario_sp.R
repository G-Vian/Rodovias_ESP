# ==============================================================================
# GRAFO RODOVIÁRIO PONDERADO DE MUNICÍPIOS PARA INLA (CAMINHO A: generic0 ICAR)
# Nós    = municípios de SP
# Arestas = rodovia ligando diretamente dois municípios (cruza a divisa)
# Peso   = quantas rodovias distintas cruzam a divisa A-B  (-> matriz W)
# Saída p/ modelo: C = D - W  (Laplaciano ponderado), usado em model="generic0"
# ------------------------------------------------------------------------------
# CAMINHOS RELATIVOS: este script localiza a si mesmo e procura as pastas de
# dados AO LADO dele (subpasta "dados_entrada/"); toda a saída vai para
# "Resultados_Grafo1/", criada automaticamente na mesma pasta do script.
# ==============================================================================

library(sf)
library(dplyr)
library(spdep)
library(Matrix)
library(ggplot2)
# library(INLA)  # necessário só na PARTE 11 (ajuste do modelo)

sf::sf_use_s2(FALSE)

# ==============================================================================
# PARTE 0: ANCORAGEM DE CAMINHOS RELATIVOS
# ------------------------------------------------------------------------------
# Descobre a pasta onde ESTE script está, de forma robusta a 3 cenários:
#   (a) Rscript caminho/script.R   (linha de comando)
#   (b) source("caminho/script.R") (console/RStudio)
#   (c) "Source" no RStudio        (botão)
# Se nada funcionar, usa o diretório de trabalho atual (getwd()).
# ==============================================================================
get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  f <- grep("^--file=", args, value = TRUE)
  if (length(f)) return(dirname(normalizePath(sub("^--file=", "", f))))
  if (!is.null(sys.frames()) && length(sys.frames())) {
    sf_path <- tryCatch(normalizePath(sys.frame(1)$ofile), error = function(e) NULL)
    if (!is.null(sf_path)) return(dirname(sf_path))
  }
  if (requireNamespace("rstudioapi", quietly = TRUE) &&
      rstudioapi::isAvailable()) {
    p <- tryCatch(rstudioapi::getActiveDocumentContext()$path, error = function(e) "")
    if (nzchar(p)) return(dirname(normalizePath(p)))
  }
  getwd()
}

DIR_SCRIPT  <- get_script_dir()
DIR_DADOS   <- file.path(DIR_SCRIPT, "dados_entrada")     # shapefiles de entrada
DIR_RESULT  <- file.path(DIR_SCRIPT, "Resultados_Grafo1") # tudo que o script gera
if (!dir.exists(DIR_RESULT)) dir.create(DIR_RESULT, recursive = TRUE)
cat("Pasta do script :", DIR_SCRIPT, "\n")
cat("Pasta de dados  :", DIR_DADOS,  "\n")
cat("Pasta de saída  :", DIR_RESULT, "\n\n")

# ==============================================================================
# PARTE 1: LEITURA  (caminhos relativos a dados_entrada/)
# ==============================================================================
path_mun   <- file.path(DIR_DADOS, "SP_Municipios_2024", "SP_Municipios_2024.shp")
path_fed   <- file.path(DIR_DADOS, "GEOFT_TRECHO_RODOVIARIO_FEDERAL", "GEOFT_TRECHO_RODOVIARIO_FEDERAL.shp")
path_est   <- file.path(DIR_DADOS, "Sistema Rodoviário Estadual", "MALHA_RODOVIARIA.shp")
path_bc250 <- file.path(DIR_DADOS, "BC250", "bc_250_shapefiles_2026_03_03", "rod_trecho_rodoviario_l.shp")

# Checagem amigável: avisa quais arquivos não foram encontrados antes de falhar
faltando <- c(path_mun, path_fed, path_est, path_bc250)[!file.exists(c(path_mun, path_fed, path_est, path_bc250))]
if (length(faltando)) {
  stop("Shapefile(s) não encontrado(s) em 'dados_entrada/':\n  ",
       paste(faltando, collapse = "\n  "),
       "\nColoque as pastas de dados ao lado do script (ver LEIA-ME).")
}

cat("Lendo shapefiles...\n")
sp_municipios <- st_read(path_mun, quiet = TRUE)
sp_malha_est  <- st_read(path_est, quiet = TRUE)
br_malha_fed  <- st_read(path_fed, quiet = TRUE)
bc250_rod     <- st_read(path_bc250, quiet = TRUE)

# ==============================================================================
# PARTE 2: CRS, VALIDAÇÃO E RECORTE
# ==============================================================================

crs_padrao <- 4674
cat("Harmonizando CRS...\n")
sp_municipios <- st_make_valid(st_transform(sp_municipios, crs_padrao))
sp_malha_est  <- st_make_valid(st_transform(sp_malha_est,  crs_padrao))
br_malha_fed  <- st_make_valid(st_transform(br_malha_fed,  crs_padrao))
bc250_rod     <- st_make_valid(st_transform(bc250_rod,     crs_padrao))

limite_sp <- st_union(sp_municipios)
cat("Recortando malha federal e BC250 para SP...\n")
sp_malha_fed <- st_collection_extract(st_intersection(br_malha_fed, limite_sp), "LINE")
bc250_sp     <- st_collection_extract(st_intersection(bc250_rod,    limite_sp), "LINE")

# ==============================================================================
# PARTE 3: REDE VIÁRIA UNIFICADA  (cada trecho recebe um id de via, rid)
# ==============================================================================

roads_all <- rbind(
  st_sf(origem = "federal",  geometry = st_geometry(sp_malha_fed)),
  st_sf(origem = "estadual", geometry = st_geometry(sp_malha_est)),
  st_sf(origem = "bc250",    geometry = st_geometry(bc250_sp))
)
roads_all <- st_make_valid(roads_all)
roads_all$rid <- seq_len(nrow(roads_all))   # identificador de via (p/ contar cruzamentos)

n_fed <- nrow(sp_malha_fed); n_est <- nrow(sp_malha_est); n_bc <- nrow(bc250_sp)

# ==============================================================================
# PARTE 4: ARESTAS + PESOS  (núcleo do Caminho A)
# ==============================================================================

sp_municipios <- sp_municipios %>% mutate(id_inla = row_number())
n <- nrow(sp_municipios)

cat("\nQuebrando rodovias nos limites municipais...\n")
pieces <- st_intersection(roads_all, sp_municipios %>% select(id_inla))
pieces <- st_collection_extract(pieces, "LINE")   # pieces carrega id_inla e rid

cat("Detectando cruzamentos de via na divisa entre municípios...\n")
rel <- st_intersects(pieces)
cross <- data.frame(i = rep(seq_along(rel), lengths(rel)), j = unlist(rel))
cross$mi <- pieces$id_inla[cross$i]
cross$mj <- pieces$id_inla[cross$j]
cross$ri <- pieces$rid[cross$i]
cross$rj <- pieces$rid[cross$j]

# Cruzamento legítimo: pedaços da MESMA via (ri==rj) em municípios DIFERENTES,
# que se tocam na divisa => essa via cruza a fronteira A-B.
cross <- cross[cross$mi != cross$mj & cross$ri == cross$rj, ]

# Peso = nº de vias DISTINTAS que cruzam a divisa de cada par (a<b)
a <- pmin(cross$mi, cross$mj)
b <- pmax(cross$mi, cross$mj)
wdf <- distinct(data.frame(a = a, b = b, rid = cross$ri))     # uma via por par
wdf <- aggregate(rid ~ a + b, data = wdf, FUN = length)       # conta vias por par
names(wdf)[3] <- "w"

cat("Pares conectados:", nrow(wdf),
    "| peso médio (nº de vias por par):", round(mean(wdf$w), 2), "\n")

# ---- TROCA DE MÉTRICA DE PESO (opcional) -------------------------------------
# Para usar EXTENSÃO de via cruzando, ou VDM/tráfego, substitua o 'w' acima:
#   - extensão: some st_length dos pedaços perto da divisa por par;
#   - tráfego : faça left_join do VDM/ARTESP pela sigla da rodovia (rid->sigla)
#               e use w = soma/máximo do VDM das vias que cruzam o par.
# A estrutura abaixo (W, C, modelo) não muda: só muda o valor de wdf$w.
# ------------------------------------------------------------------------------

# Matriz de pesos W (simétrica, esparsa)
W <- sparseMatrix(i = wdf$a, j = wdf$b, x = wdf$w, dims = c(n, n))
W <- W + t(W)                       # simetriza (a<b, então sem diagonal)

# Lista de vizinhos binária (p/ .graph, diagnóstico e visualização)
Abin  <- as(W > 0, "dMatrix")
nb_sp <- mat2listw(Abin, style = "B", zero.policy = TRUE)$neighbours
attr(nb_sp, "region.id") <- as.character(sp_municipios$id_inla)

# Estado ANTES do tratamento de isolados (para o log)
graus_pre   <- card(nb_sp)
comp_pre    <- n.comp.nb(nb_sp)$nc
iso         <- which(graus_pre == 0)
nomes_iso   <- sp_municipios$NM_MUN[iso]

# ==============================================================================
# PARTE 5: TRATAMENTO DE ISOLADOS + DIAGNÓSTICO DE GRAUS
# ==============================================================================

cat("\nNós isolados (grau 0):", length(iso),
    "| componentes:", comp_pre, "\n")

arestas_artificiais <- data.frame(de = character(0), para = character(0))
if (length(iso) > 0) {
  cat("Conectando", length(iso), "isolado(s) ao município mais próximo (peso 1)...\n")
  cent_m <- st_centroid(st_geometry(st_transform(sp_municipios, 31983)))
  for (k in iso) {
    d  <- st_distance(cent_m[k], cent_m[-k])
    nn <- (seq_len(n))[-k][which.min(d)]
    nb_sp[[k]]  <- sort(unique(c(setdiff(nb_sp[[k]],  0L), as.integer(nn))))
    nb_sp[[nn]] <- sort(unique(c(setdiff(nb_sp[[nn]], 0L), as.integer(k))))
    W[k, nn] <- 1; W[nn, k] <- 1     # aresta artificial também entra em W
    arestas_artificiais <- rbind(arestas_artificiais,
                                 data.frame(de = sp_municipios$NM_MUN[k], para = sp_municipios$NM_MUN[nn]))
    cat("  -", sp_municipios$NM_MUN[k], "->", sp_municipios$NM_MUN[nn], "\n")
  }
}
stopifnot(sum(card(nb_sp) == 0) == 0)

comp   <- n.comp.nb(nb_sp)
ncomp  <- comp$nc
compid <- comp$comp.id

# ---- Diagnóstico de conectividade --------------------------------------------
graus  <- card(nb_sp)                 # grau (nº de vizinhos) após correções
forca  <- rowSums(W)                  # grau ponderado (soma dos pesos das arestas)

# Municípios pouco conectados (candidatos a sub-conexão urbana / escala BC250)
idx_baixo   <- which(graus <= 2)
nomes_baixo <- sp_municipios$NM_MUN[idx_baixo]

# Cruzamento com a Região Metropolitana de SP (39 municípios) -> sinal de alerta
rmsp <- c("Arujá","Barueri","Biritiba Mirim","Caieiras","Cajamar","Carapicuíba",
          "Cotia","Diadema","Embu das Artes","Embu-Guaçu","Ferraz de Vasconcelos",
          "Francisco Morato","Franco da Rocha","Guararema","Guarulhos",
          "Itapecerica da Serra","Itapevi","Itaquaquecetuba","Jandira","Juquitiba",
          "Mairiporã","Mauá","Mogi das Cruzes","Osasco","Pirapora do Bom Jesus",
          "Poá","Ribeirão Pires","Rio Grande da Serra","Salesópolis","Santa Isabel",
          "Santana de Parnaíba","Santo André","São Bernardo do Campo",
          "São Caetano do Sul","São Lourenço da Serra","São Paulo","Suzano",
          "Taboão da Serra","Vargem Grande Paulista")
baixo_rmsp <- intersect(nomes_baixo, rmsp)

# Mais conectados (por grau e por força)
ord_grau  <- order(graus, decreasing = TRUE)
ord_forca <- order(forca, decreasing = TRUE)
top_grau  <- data.frame(mun = sp_municipios$NM_MUN[head(ord_grau, 10)],
                        grau = graus[head(ord_grau, 10)])
top_forca <- data.frame(mun = sp_municipios$NM_MUN[head(ord_forca, 10)],
                        forca = forca[head(ord_forca, 10)])

cat("\n---- DISTRIBUIÇÃO DE GRAUS ----\n"); print(summary(graus))
cat("Municípios com grau <= 2:", length(idx_baixo), "\n")
if (length(baixo_rmsp) > 0)
  cat("DESTES, pertencem à RMSP (alerta de sub-conexão urbana):",
      paste(baixo_rmsp, collapse = ", "), "\n")
cat("Componentes conexas finais:", ncomp, "\n")

# ==============================================================================
# PARTE 6: MATRIZ DE ESTRUTURA  C = D - W  (Laplaciano ponderado -> generic0)
# ==============================================================================

D <- Diagonal(x = rowSums(W))
C <- D - W                          # simétrica, PSD, singular (rankdef = ncomp)
C <- as(C, "TsparseMatrix")

# Restrições soma-zero: UMA por componente conexa (ICAR exige isso)
Aconstr <- matrix(0, nrow = ncomp, ncol = n)
for (cc in seq_len(ncomp)) Aconstr[cc, compid == cc] <- 1
econstr <- rep(0, ncomp)

# ==============================================================================
# PARTE 7: EXPORTAÇÃO (grafo binário + estrutura ponderada)  -> Resultados_Grafo1/
# ==============================================================================

nb2INLA(file.path(DIR_RESULT, "sp_grafo_rodoviario.graph"), nb_sp)   # binário (referência)
saveRDS(W,       file.path(DIR_RESULT, "W_pesos.rds"))               # matriz de pesos
saveRDS(C,       file.path(DIR_RESULT, "C_estrutura.rds"))           # matriz p/ generic0
saveRDS(Aconstr, file.path(DIR_RESULT, "Aconstr.rds"))
saveRDS(econstr, file.path(DIR_RESULT, "econstr.rds"))
cat("\nW, C e restrições salvas em:", DIR_RESULT, "\n")

# ==============================================================================
# PARTE 8: MAPA DO GRAFO (espessura da aresta ~ peso)
# ==============================================================================

cent   <- st_centroid(st_geometry(sp_municipios))
coords <- st_coordinates(cent)
Wt     <- as(triu(W), "TsparseMatrix")          # triângulo superior p/ não duplicar
edge_df <- data.frame(from = Wt@i + 1, to = Wt@j + 1, w = Wt@x)
edge_df <- edge_df[edge_df$from != edge_df$to, ]

edge_lines <- st_sf(
  w = edge_df$w,
  geometry = st_sfc(lapply(seq_len(nrow(edge_df)), function(r)
    st_linestring(rbind(coords[edge_df$from[r], ], coords[edge_df$to[r], ]))),
    crs = crs_padrao)
)

mapa_grafo <- ggplot() +
  geom_sf(data = sp_municipios, fill = "gray97", color = "gray80", size = 0.1) +
  geom_sf(data = edge_lines, aes(linewidth = w), color = "steelblue", alpha = 0.6) +
  scale_linewidth(range = c(0.2, 1.6), name = "nº de vias") +
  geom_sf(data = cent, color = "firebrick", size = 0.5) +
  theme_minimal() +
  labs(title = "Grafo rodoviário ponderado - SP",
       subtitle = "Espessura da aresta = nº de rodovias ligando os municípios",
       caption = "Fontes: BC250 (IBGE 2025), DER-SP, DNIT")
print(mapa_grafo)
ggsave(file.path(DIR_RESULT, "Grafo_Rodoviario_Ponderado_SP.png"),
       mapa_grafo, width = 12, height = 8, dpi = 300, bg = "white")

# ==============================================================================
# PARTE 8B: MAPA CARTOGRÁFICO DAS RODOVIAS REAIS (traçado, por fonte)
# ------------------------------------------------------------------------------
# Diferente do grafo (linhas retas centroide-a-centroide), este mostra a
# GEOMETRIA REAL das vias — útil p/ conferir a cobertura de cada fonte viária.
# ==============================================================================
mapa_rodovias <- ggplot() +
  geom_sf(data = sp_municipios, fill = "gray98", color = "gray85", linewidth = 0.1) +
  geom_sf(data = roads_all, aes(color = origem), linewidth = 0.3, alpha = 0.8) +
  scale_color_manual(values = c(federal = "firebrick", estadual = "steelblue", bc250 = "darkgreen"),
                     name = "fonte da via") +
  theme_minimal() +
  labs(title = "Malha rodoviária de SP (traçado real)",
       subtitle = "Federal (DNIT) + Estadual (DER) + BC250 (IBGE)",
       caption = "Geometria real das vias usadas para construir o grafo")
ggsave(file.path(DIR_RESULT, "Mapa_Rodovias_SP.png"),
       mapa_rodovias, width = 12, height = 8, dpi = 300, bg = "white")

# ==============================================================================
# PARTE 9: LOG DETALHADO  -> Resultados_Grafo1/
# ==============================================================================

linha <- function() cat(strrep("-", 72), "\n")

sink(file.path(DIR_RESULT, "LOG_Grafo_Ponderado.txt"))

cat("########################################################################\n")
cat("RELATÓRIO - GRAFO RODOVIÁRIO PONDERADO DE MUNICÍPIOS (SP)\n")
cat("Projeto: Epidemiologia Matemática / INLA\n")
cat("Gerado em:", as.character(Sys.time()), "\n")
cat("########################################################################\n\n")

cat("[1] DEFINIÇÃO DO GRAFO\n"); linha()
cat("- Nó    : município (", n, "no total; ordem = id_inla = linha de sp_municipios).\n")
cat("- Aresta: existe rodovia (federal/estadual/BC250) cruzando a divisa A-B.\n")
cat("- Peso  : nº de rodovias DISTINTAS que cruzam a divisa do par.\n")
cat("- CRS de trabalho (topologia): EPSG", crs_padrao, "(SIRGAS2000 geográfico).\n")
cat("- Saída p/ modelo: C = D - W (Laplaciano ponderado), p/ model='generic0'.\n\n")

cat("[2] INSUMOS VIÁRIOS (após recorte para SP)\n"); linha()
cat("- Trechos federais (DNIT) :", n_fed, "\n")
cat("- Trechos estaduais (DER) :", n_est, "\n")
cat("- Trechos BC250 (IBGE)    :", n_bc, "\n")
cat("- Total de trechos (rid)  :", nrow(roads_all), "\n")
cat("- Pedaços após recorte municipal:", nrow(pieces), "\n\n")

cat("[3] ESTATÍSTICAS DAS ARESTAS E PESOS\n"); linha()
cat("- Pares de municípios conectados :", nrow(wdf), "\n")
cat("- Peso médio (vias por par)      :", round(mean(wdf$w), 2), "\n")
cat("- Peso mediano                   :", median(wdf$w), "\n")
cat("- Peso máximo                    :", max(wdf$w), "\n")
cat("- Peso mínimo                    :", min(wdf$w), "\n")
cat("- Desvio-padrão do peso          :", round(sd(wdf$w), 2), "\n")
cat("- Distribuição do peso (quantis):\n")
print(quantile(wdf$w, probs = c(0,.25,.5,.75,.9,.95,1)))
cat("\n")

cat("[4] CONECTIVIDADE DOS NÓS (após correções)\n"); linha()
cat("- Distribuição do GRAU (nº de vizinhos):\n"); print(summary(graus))
cat("- Distribuição da FORÇA (soma de pesos):\n"); print(summary(forca))
cat("\nTop 10 municípios mais conectados (por grau):\n")
print(top_grau, row.names = FALSE)
cat("\nTop 10 municípios mais conectados (por força = soma de pesos):\n")
print(top_forca, row.names = FALSE)
cat("\n")

cat("[5] ALERTA: SUB-CONEXÃO EM ÁREA URBANA (escala da BC250)\n"); linha()
cat("A BC250 é escala 1:250.000 e NAO representa malha de ruas urbanas. Em áreas\n")
cat("densamente urbanizadas (RMSP), vias de passagem podem ter sido atribuídas aos\n")
cat("vizinhos e nenhuma rodovia federal/estadual cruza a divisa de municípios\n")
cat("pequenos. Resultado: municípios urbanos pequenos podem ficar SUB-CONECTADOS\n")
cat("(ou até isolados) no grafo, justamente onde a conectividade real é mais densa.\n\n")
cat("Exemplo observado: Ferraz de Vasconcelos (município urbano da Grande SP,\n")
cat("cercado de vizinhos) saiu como grau 0 no grafo puramente rodoviário.\n\n")
cat("- Municípios com grau <= 2:", length(idx_baixo), "\n")
if (length(idx_baixo) > 0) {
  cat("  Lista:\n")
  cat("  ", paste(nomes_baixo, collapse = ", "), "\n")
}
cat("- DESTES, pertencentes à RMSP (alerta forte):", length(baixo_rmsp), "\n")
if (length(baixo_rmsp) > 0) {
  cat("  ", paste(baixo_rmsp, collapse = ", "), "\n")
  cat("\n  >> Se há muitos nomes da Grande SP acima, o grafo rodoviário está\n")
  cat("     subestimando a conectividade metropolitana. Considere um grafo\n")
  cat("     HÍBRIDO (rodoviário ponderado + contiguidade poly2nb onde faltar).\n")
}
cat("\n")

cat("[6] NÓS ISOLADOS E ARESTAS ARTIFICIAIS\n"); linha()
cat("- Componentes conexas ANTES da correção:", comp_pre, "\n")
cat("- Nós isolados (grau 0) ANTES da correção:", length(iso), "\n")
if (length(iso) > 0) cat("  ", paste(nomes_iso, collapse = ", "), "\n")
cat("- Arestas artificiais criadas (ligação ao vizinho mais próximo, peso 1):\n")
if (nrow(arestas_artificiais) > 0) {
  for (r in seq_len(nrow(arestas_artificiais)))
    cat("    *", arestas_artificiais$de[r], "->", arestas_artificiais$para[r], "\n")
  cat("\n  NOTA: a aresta de Ilhabela foi criada por proximidade de centroide.\n")
  cat("  Se quiser fidelidade à balsa real (Ilhabela <-> São Sebastião), force-a\n")
  cat("  manualmente em W antes de montar C.\n")
}
cat("- Componentes conexas APÓS a correção:", ncomp, "\n\n")

cat("[7] ARTEFATOS GERADOS (em Resultados_Grafo1/)\n"); linha()
cat("- sp_grafo_rodoviario.graph  (grafo binário, referência/diagnóstico)\n")
cat("- W_pesos.rds                (matriz de pesos simétrica)\n")
cat("- C_estrutura.rds            (C = D - W, para model='generic0')\n")
cat("- Aconstr.rds / econstr.rds  (restrições soma-zero, 1 por componente)\n")
cat("- Grafo_Rodoviario_Ponderado_SP.png (diagrama do grafo)\n")
cat("- Mapa_Rodovias_SP.png       (traçado real das rodovias por fonte)\n\n")

cat("[8] PRÓXIMO PASSO (INLA)\n"); linha()
cat("Caminho A: f(id_inla, model='generic0', Cmatrix=C_scaled,\n")
cat("            constr=FALSE, extraconstr=list(A=Aconstr,e=econstr),\n")
cat("            rankdef=", ncomp, ")\n", sep = "")
cat("Verossimilhança: family='nbinomial'. 'dados' deve estar ordenado por id_inla.\n")

sink()
cat("Concluído. Log detalhado salvo em:", file.path(DIR_RESULT, "LOG_Grafo_Ponderado.txt"), "\n")

# ==============================================================================
# PARTE 10/11: AJUSTE DO MODELO NO INLA  (Caminho A)  -- descomente p/ rodar
# ==============================================================================
# library(INLA)
#
# # 'dados' deve ter UMA linha por município, na MESMA ORDEM de sp_municipios,
# # com: id_inla, casos (contagem) e E (esperado/offset, ex.: pop * taxa base).
# # Ex.: dados <- sp_municipios %>% st_drop_geometry() %>%
# #               left_join(meus_casos, by = "CD_MUN") %>% arrange(id_inla)
#
# C_scaled <- inla.scale.model(C, constr = list(A = Aconstr, e = econstr))
#
# formula <- casos ~ 1 +
#   f(id_inla,
#     model       = "generic0",
#     Cmatrix     = C_scaled,
#     constr      = FALSE,
#     extraconstr = list(A = Aconstr, e = econstr),
#     rankdef     = ncomp,
#     hyper       = list(prec = list(prior = "pc.prec", param = c(1, 0.01))))
#
# fit <- inla(formula,
#             family            = "nbinomial",
#             data              = dados,
#             E                 = dados$E,
#             control.predictor = list(compute = TRUE),
#             control.compute   = list(dic = TRUE, waic = TRUE))
# summary(fit)