# ==============================================================================
# GRAFO #2 - GRAFO RODOVIÁRIO COM PESOS FUNDIDOS (4 FONTES)
# Mesma TOPOLOGIA do Grafo #1; peso de aresta = fusão de:
#   (a) nº de vias  (b) VDM/DER  (c) OD/Metrô  (d) IBGE VAR07
# Normalização: log1p + z-score | Combinação: média ponderada das fontes PRESENTES
# Saída: C = D - W_fundida (para model="generic0" no INLA, family="nbinomial")
# ------------------------------------------------------------------------------
# CAMINHOS RELATIVOS: o script se localiza sozinho. Ele LÊ insumos gerados por
# outros scripts (W_pesos.rds do Grafo #1; CSVs do inventário) e a saída vai
# para "Resultados_Grafo2/". Veja o bloco de CONFIGURAÇÃO abaixo.
# Pré-requisito: rodar antes o Grafo #1 (gera W_pesos.rds) e o inventário (CSVs).
# ==============================================================================

library(sf)
library(dplyr)
library(Matrix)
library(spdep)
library(ggplot2)

sf::sf_use_s2(FALSE)

# ==============================================================================
# PARTE 0: ANCORAGEM DE CAMINHOS RELATIVOS
# ==============================================================================
get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  f <- grep("^--file=", args, value = TRUE)
  if (length(f)) return(dirname(normalizePath(sub("^--file=", "", f))))
  if (!is.null(sys.frames()) && length(sys.frames())) {
    sf_path <- tryCatch(normalizePath(sys.frame(1)$ofile), error = function(e) NULL)
    if (!is.null(sf_path)) return(dirname(sf_path))
  }
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    p <- tryCatch(rstudioapi::getActiveDocumentContext()$path, error = function(e) "")
    if (nzchar(p)) return(dirname(normalizePath(p)))
  }
  getwd()
}
DIR_SCRIPT <- get_script_dir()
DIR_RESULT <- file.path(DIR_SCRIPT, "Resultados_Grafo2")
if (!dir.exists(DIR_RESULT)) dir.create(DIR_RESULT, recursive = TRUE)

# ------------------------------------------------------------------------------
# >>> CONFIGURAÇÃO (mude aqui) <<<
# ------------------------------------------------------------------------------
pesos_fonte <- c(vias = 1, vdm = 1, od = 1, ibge = 1)  # importância relativa
metodo_norm <- "zscore"   # "zscore" (log1p+z) ou "rank" (posto 0-1)
floor_peso  <- 0.05       # piso do peso final (>0, evita aresta com peso 0)

# --- ONDE ESTÃO OS INSUMOS (cada um mora num lugar diferente) ---
# Dados ORIGINAIS ficam em dados_entrada/.
# SAÍDAS de scripts anteriores ficam nas pastas de resultados deles.
DIR_DADOS  <- file.path(DIR_SCRIPT, "dados_entrada")        # shapefiles originais
DIR_W_VIAS <- file.path(DIR_SCRIPT, "Resultados_Grafo1")    # W_pesos.rds (Grafo #1)
DIR_CSVS   <- file.path(DIR_SCRIPT, "Resultados_Inventario")# CSVs (inventário)
PASTA_OD   <- "Origem_Destinos_2023"                        # dicionário OD (em dados_entrada/)

# --- Caminhos montados a partir do acima ---
path_mun     <- file.path(DIR_DADOS, "SP_Municipios_2024", "SP_Municipios_2024.shp")
path_W_vias  <- file.path(DIR_W_VIAS, "W_pesos.rds")          # do Grafo #1
path_vdm     <- file.path(DIR_CSVS,   "der_vdm_pares.csv")
path_od      <- file.path(DIR_CSVS,   "od_municipios_fluxo.csv")
path_ibge    <- file.path(DIR_CSVS,   "ibge_ligacoes_sp.csv")
path_od_mun  <- file.path(DIR_DADOS, PASTA_OD, "Site_190225_PesquisaOD2023", "Site_190225",
                          "002_Site Metro Mapas_190225", "Shape", "Municipios_2023.shp")

# Checagem amigável dos insumos obrigatórios
.req <- c(municipios = path_mun, W_vias = path_W_vias,
          vdm = path_vdm, od = path_od, ibge = path_ibge)
.falta <- .req[!file.exists(.req)]
if (length(.falta)) {
  stop("Insumo(s) não encontrado(s):\n  ",
       paste(sprintf("[%s] %s", names(.falta), .falta), collapse = "\n  "),
       "\nAjuste DIR_W_VIAS / DIR_CSVS / nomes de pasta no bloco de CONFIGURAÇÃO.")
}

# ==============================================================================
# PARTE 1: BASE - municípios, id_inla, e arestas do Grafo #1
# ==============================================================================
sp_municipios <- st_read(path_mun, quiet = TRUE)
sp_municipios <- sp_municipios %>% mutate(id_inla = row_number())
n  <- nrow(sp_municipios)
cd <- as.character(sp_municipios$CD_MUN)   # código IBGE 7 dígitos
nm <- as.character(sp_municipios$NM_MUN)

W_vias <- readRDS(path_W_vias)             # matriz simétrica (nº de vias) do Grafo #1
stopifnot(nrow(W_vias) == n)               # ordem TEM que bater com sp_municipios

Wt <- as(Matrix::triu(W_vias), "TsparseMatrix")
ed <- data.frame(iA = Wt@i + 1L, iB = Wt@j + 1L, vias = Wt@x)
ed <- ed[ed$iA != ed$iB & ed$vias > 0, ]
ed$cdA <- cd[ed$iA]; ed$cdB <- cd[ed$iB]
ed$key <- paste(pmin(ed$cdA, ed$cdB), pmax(ed$cdA, ed$cdB), sep = "_")
cat("Arestas do Grafo #1:", nrow(ed), "\n")

mk_key <- function(a, b) paste(pmin(a, b), pmax(a, b), sep = "_")
nrm <- function(x) {
  x <- toupper(trimws(x)); x <- iconv(x, to = "ASCII//TRANSLIT")
  x <- gsub("[^A-Z ]", "", x); trimws(gsub("\\s+", " ", x))
}

# ==============================================================================
# PARTE 2: FONTE VDM (DER) - por nome de município
# ==============================================================================
ed$vdm <- NA_real_
tryCatch({
  vdm <- read.csv(path_vdm, stringsAsFactors = FALSE, fileEncoding = "UTF-8")
  vdm <- vdm[vdm$split_ok %in% c(TRUE, "TRUE") & !is.na(vdm$VDM2024), ]
  nm_norm <- nrm(nm)
  dup <- unique(nm_norm[duplicated(nm_norm)])
  if (length(dup) > 0) cat("AVISO: nomes normalizados duplicados:", paste(dup, collapse=", "), "\n")
  lut <- setNames(cd, nm_norm)
  vdm$cdA <- lut[nrm(vdm$mun_A)]; vdm$cdB <- lut[nrm(vdm$mun_B)]
  nao_casou <- vdm[is.na(vdm$cdA) | is.na(vdm$cdB), c("mun_A","mun_B")]
  vdm <- vdm[!is.na(vdm$cdA) & !is.na(vdm$cdB), ]
  vdm$key <- mk_key(vdm$cdA, vdm$cdB)
  vdm_k <- aggregate(VDM2024 ~ key, vdm, sum)
  ed$vdm <- vdm_k$VDM2024[match(ed$key, vdm_k$key)]
  cat("VDM: pares casados a arestas:", sum(!is.na(ed$vdm)),
      "| pares VDM sem nome casado:", nrow(nao_casou), "\n")
  if (nrow(nao_casou) > 0) { cat("  não casados (revisar):\n"); print(head(nao_casou, 20)) }
}, error = function(e) cat("ERRO VDM:", conditionMessage(e), "\n"))

# ==============================================================================
# PARTE 3: FONTE OD (Metrô) - código OD -> CD_IBGE pelo dicionário
# ==============================================================================
ed$od <- NA_real_
tryCatch({
  mun_od <- st_read(path_od_mun, quiet = TRUE)
  dict <- setNames(as.character(mun_od$CD_IBGE), as.character(mun_od$NumeroMuni))
  od <- read.csv(path_od, stringsAsFactors = FALSE)
  od$cdO <- dict[as.character(od$MUNI_O)]; od$cdD <- dict[as.character(od$MUNI_D)]
  fora <- sum(is.na(od$cdO) | is.na(od$cdD))
  od <- od[!is.na(od$cdO) & !is.na(od$cdD), ]
  od$key <- mk_key(od$cdO, od$cdD)
  od_k <- aggregate(viagens ~ key, od, sum)
  ed$od <- od_k$viagens[match(ed$key, od_k$key)]
  cat("OD: pares casados a arestas:", sum(!is.na(ed$od)),
      "| registros OD sem código no dicionário:", fora, "\n")
}, error = function(e) cat("ERRO OD:", conditionMessage(e), "\n"))

# ==============================================================================
# PARTE 4: FONTE IBGE (VAR07 = frequência total de saídas) - join por CD direto
# ==============================================================================
ed$ibge <- NA_real_
tryCatch({
  ib <- read.csv(path_ibge, stringsAsFactors = FALSE, fileEncoding = "UTF-8")
  ib$cdA <- as.character(ib$CODMUNDV_A); ib$cdB <- as.character(ib$CODMUNDV_B)
  ib$key <- mk_key(ib$cdA, ib$cdB)
  ib_k <- aggregate(VAR07 ~ key, ib, sum)
  ed$ibge <- ib_k$VAR07[match(ed$key, ib_k$key)]
  cat("IBGE: pares casados a arestas:", sum(!is.na(ed$ibge)), "\n")
}, error = function(e) cat("ERRO IBGE:", conditionMessage(e), "\n"))

# ==============================================================================
# PARTE 5: NORMALIZAÇÃO + FUSÃO + POSITIVAÇÃO
# ==============================================================================
normaliza <- function(x, metodo) {
  ok <- !is.na(x); out <- rep(NA_real_, length(x))
  if (sum(ok) < 2) return(out)
  if (metodo == "zscore") {
    lx <- log1p(x[ok]); s <- sd(lx)
    out[ok] <- if (s > 0) (lx - mean(lx)) / s else 0
  } else { out[ok] <- (rank(x[ok]) - 1) / (sum(ok) - 1) }
  out
}
Z <- cbind(vias = normaliza(ed$vias, metodo_norm), vdm = normaliza(ed$vdm, metodo_norm),
           od = normaliza(ed$od, metodo_norm),     ibge = normaliza(ed$ibge, metodo_norm))
pw <- pesos_fonte[colnames(Z)]
Wt_mat <- matrix(pw, nrow(Z), ncol(Z), byrow = TRUE); Wt_mat[is.na(Z)] <- NA
fundido <- rowSums(Z * Wt_mat, na.rm = TRUE) / rowSums(Wt_mat, na.rm = TRUE)
ed$n_fontes <- rowSums(!is.na(Z))
rng <- range(fundido, na.rm = TRUE)
ed$w_final <- if (diff(rng) > 0)
  floor_peso + (1 - floor_peso) * (fundido - rng[1]) / diff(rng) else 1

# ==============================================================================
# PARTE 6: MATRIZ W_fundida, C = D - W, restrições  (para generic0)
# ==============================================================================
Wf <- sparseMatrix(i = ed$iA, j = ed$iB, x = ed$w_final, dims = c(n, n)); Wf <- Wf + t(Wf)
nb_chk <- mat2listw(as(Wf > 0, "dMatrix"), style = "B", zero.policy = TRUE)$neighbours
comp <- n.comp.nb(nb_chk); ncomp <- comp$nc; compid <- comp$comp.id
stopifnot(sum(card(nb_chk) == 0) == 0)
D <- Diagonal(x = rowSums(Wf)); C <- as(D - Wf, "TsparseMatrix")
Aconstr <- matrix(0, nrow = ncomp, ncol = n)
for (cc in seq_len(ncomp)) Aconstr[cc, compid == cc] <- 1
econstr <- rep(0, ncomp)

# ==============================================================================
# PARTE 7: EXPORTAÇÃO (matrizes + tabela de arestas)  -> Resultados_Grafo2/
# ==============================================================================
saveRDS(Wf,      file.path(DIR_RESULT, "W_fundida.rds"))
saveRDS(C,       file.path(DIR_RESULT, "C_fundida.rds"))
saveRDS(Aconstr, file.path(DIR_RESULT, "Aconstr_fundida.rds"))
saveRDS(econstr, file.path(DIR_RESULT, "econstr_fundida.rds"))
nb2INLA(file.path(DIR_RESULT, "sp_grafo2_fundido.graph"), nb_chk)
tab <- data.frame(
  munA = nm[ed$iA], munB = nm[ed$iB], cdA = ed$cdA, cdB = ed$cdB,
  vias = ed$vias, vdm = ed$vdm, od = ed$od, ibge = ed$ibge,
  z_vias = Z[, "vias"], z_vdm = Z[, "vdm"], z_od = Z[, "od"], z_ibge = Z[, "ibge"],
  n_fontes = ed$n_fontes, peso_final = ed$w_final)
write.csv(tab, file.path(DIR_RESULT, "grafo2_arestas_pesos.csv"), row.names = FALSE, fileEncoding = "UTF-8")

# ==============================================================================
# PARTE 8: VISUALIZAÇÃO DO GRAFO  (2 mapas em PNG)
#  Mapa A: largura/cor da aresta ~ peso fundido (intensidade da ligação)
#  Mapa B: cor da aresta ~ nº de fontes (revela onde há fluxo VDM/OD)
# ==============================================================================
cent   <- suppressWarnings(st_centroid(st_geometry(sp_municipios)))
coords <- st_coordinates(cent)

edge_sf <- st_sf(
  peso     = ed$w_final,
  n_fontes = factor(ed$n_fontes, levels = sort(unique(ed$n_fontes))),
  geometry = st_sfc(lapply(seq_len(nrow(ed)), function(r)
    st_linestring(rbind(coords[ed$iA[r], ], coords[ed$iB[r], ]))),
    crs = st_crs(sp_municipios))
)
edge_sf <- edge_sf[order(edge_sf$peso), ]   # arestas fortes desenhadas por cima

# ---- Mapa A: peso fundido ----
mapa_peso <- ggplot() +
  geom_sf(data = sp_municipios, fill = "gray98", color = "gray85", linewidth = 0.1) +
  geom_sf(data = edge_sf, aes(linewidth = peso, color = peso), alpha = 0.85) +
  scale_linewidth(range = c(0.15, 1.8), guide = "none") +
  scale_color_viridis_c(option = "C", name = "peso\nfundido") +
  geom_sf(data = cent, color = "black", size = 0.25) +
  theme_minimal() +
  labs(title = "Grafo #2 - peso fundido das arestas",
       subtitle = "Cor/largura = intensidade combinada (vias + VDM + OD + IBGE)",
       caption = "Fontes: BC250/DER/DNIT (vias), DER (VDM), Metrô (OD), IBGE (VAR07)")
png(file.path(DIR_RESULT, "Grafo2_mapa_peso.png"), width = 12, height = 9, units = "in", res = 300, bg = "white")
print(mapa_peso)
dev.off()
# ---- Mapa B: cobertura de fontes (mostra concentração do fluxo na RMSP) ----
cores_fontes <- c("1" = "grey75", "2" = "steelblue", "3" = "firebrick", "4" = "darkgreen")
mapa_fontes <- ggplot() +
  geom_sf(data = sp_municipios, fill = "gray98", color = "gray85", linewidth = 0.1) +
  geom_sf(data = edge_sf, aes(color = n_fontes), linewidth = 0.5, alpha = 0.85) +
  scale_color_manual(values = cores_fontes, name = "nº de fontes\nna aresta") +
  theme_minimal() +
  labs(title = "Grafo #2 - nº de fontes por aresta",
       subtitle = "Azul/vermelho = arestas com VDM/OD (concentradas na Grande SP)",
       caption = "Cinza = só nº de vias | a maioria das arestas tem 1-2 fontes")
png(file.path(DIR_RESULT, "Grafo2_mapa_fontes.png"), width = 12, height = 9, units = "in", res = 300, bg = "white")
print(mapa_fontes)
dev.off()
print(mapa_peso); print(mapa_fontes)
cat("Mapas salvos: Grafo2_mapa_peso.png e Grafo2_mapa_fontes.png\n")

# ==============================================================================
# PARTE 9: LOG / RELATÓRIO DE EXECUÇÃO  -> Resultados_Grafo2/
# ==============================================================================
sink(file.path(DIR_RESULT, "LOG_Grafo2_Fusao.txt"))
cat("GRAFO #2 - PESOS FUNDIDOS - LOG DE EXECUÇÃO\n")
cat("Gerado em:", as.character(Sys.time()), "\n\n")
cat("Config: metodo_norm =", metodo_norm, "| floor_peso =", floor_peso, "\n")
cat("Pesos por fonte:", paste(names(pesos_fonte), pesos_fonte, sep="=", collapse=" "), "\n\n")
cat("Arestas (topologia do Grafo #1):", nrow(ed), "\n")
cat("Cobertura por fonte (nº de arestas com a fonte):\n")
cat("  vias :", sum(!is.na(ed$vias)), "(sempre presente)\n")
cat("  VDM  :", sum(!is.na(ed$vdm)),  "\n")
cat("  OD   :", sum(!is.na(ed$od)),   "\n")
cat("  IBGE :", sum(!is.na(ed$ibge)), "\n\n")
cat("Distribuição do nº de fontes por aresta:\n"); print(table(ed$n_fontes))
cat("\nResumo do peso final:\n"); print(summary(ed$w_final))
cat("\nComponentes conexas:", ncomp, "(deve ser 1, igual ao Grafo #1)\n")
cat("\nArtefatos salvos (em Resultados_Grafo2/):\n")
cat("  W_fundida.rds, C_fundida.rds, Aconstr_fundida.rds, econstr_fundida.rds\n")
cat("  sp_grafo2_fundido.graph, grafo2_arestas_pesos.csv\n")
cat("  Grafo2_mapa_peso.png, Grafo2_mapa_fontes.png\n")
sink()

cat("\nGrafo #2 concluído.\n")
cat("Cobertura -> VDM:", sum(!is.na(ed$vdm)), "| OD:", sum(!is.na(ed$od)),
    "| IBGE:", sum(!is.na(ed$ibge)), "| arestas:", nrow(ed), "\n")

# ==============================================================================
# PARTE 10: USO NO INLA  (Caminho A) -- descomente para rodar
# ==============================================================================
# library(INLA)
# C_scaled <- inla.scale.model(C, constr = list(A = Aconstr, e = econstr))
# formula <- casos ~ 1 +
#   f(id_inla, model="generic0", Cmatrix=C_scaled, constr=FALSE,
#     extraconstr=list(A=Aconstr, e=econstr), rankdef=ncomp,
#     hyper=list(prec=list(prior="pc.prec", param=c(1,0.01))))
# fit <- inla(formula, family="nbinomial", data=dados, E=dados$E,
#             control.compute=list(dic=TRUE, waic=TRUE))