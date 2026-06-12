# ==============================================================================
# GRAFO #3 - ADJACÊNCIA + ATALHOS  (RÉGUA ÚNICA DE NORMALIZAÇÃO)
# Topologia: base (vizinhança rodoviária) + atalhos de proximidade + saltos OD.
# DIFERENÇA-CHAVE p/ versão anterior: TODAS as arestas e fontes são normalizadas
# num ÚNICO cálculo (uma régua por fonte, global), garantindo comparabilidade.
# Por isso lê os valores CRUS (W_pesos.rds = vias) e RE-JUNTA VDM/OD/IBGE.
# Saída: C = D - W3 (para model="generic0" no INLA, family="nbinomial").
# ------------------------------------------------------------------------------
# CAMINHOS RELATIVOS: o script se localiza sozinho. LÊ insumos de pastas
# organizadas (dados_entrada/, Resultados_Grafo1/, Resultados_Inventario/) e a
# saída vai para "Resultados_Grafo3/". Veja o bloco de CONFIGURAÇÃO.
# Pré-requisitos: Grafo #1 (gera W_pesos.rds) e inventário (CSVs de VDM/OD/IBGE).
# ==============================================================================

library(sf); library(dplyr); library(Matrix); library(spdep); library(ggplot2)
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
DIR_RESULT <- file.path(DIR_SCRIPT, "Resultados_Grafo3")
if (!dir.exists(DIR_RESULT)) dir.create(DIR_RESULT, recursive = TRUE)

# ------------------------------------------------------------------------------
# >>> CONFIGURAÇÃO <<<
# ------------------------------------------------------------------------------
TETO_KM     <- 150
K_VIZ       <- 8
floor_peso  <- 0.05
pesos_fonte <- c(vias = 1, vdm = 1, od = 1, ibge = 1, prox = 1)  # importância por fonte

# --- ONDE ESTÃO OS INSUMOS (cada um na sua pasta) ---
DIR_DADOS  <- file.path(DIR_SCRIPT, "dados_entrada")         # shapefiles originais
DIR_W_VIAS <- file.path(DIR_SCRIPT, "Resultados_Grafo1")     # W_pesos.rds (Grafo #1)
DIR_CSVS   <- file.path(DIR_SCRIPT, "Resultados_Inventario") # CSVs (inventário)
PASTA_OD   <- "Origem_Destinos_2023"                         # dicionário OD (em dados_entrada/)

path_mun    <- file.path(DIR_DADOS, "SP_Municipios_2024", "SP_Municipios_2024.shp")
path_Wvias  <- file.path(DIR_W_VIAS, "W_pesos.rds")          # VIAS CRUAS (Grafo #1)
path_vdm    <- file.path(DIR_CSVS,   "der_vdm_pares.csv")
path_od     <- file.path(DIR_CSVS,   "od_municipios_fluxo.csv")
path_ibge   <- file.path(DIR_CSVS,   "ibge_ligacoes_sp.csv")
path_od_mun <- file.path(DIR_DADOS, PASTA_OD, "Site_190225_PesquisaOD2023", "Site_190225",
                         "002_Site Metro Mapas_190225", "Shape", "Municipios_2023.shp")

# Checagem amigável dos insumos obrigatórios
.req <- c(municipios = path_mun, W_vias = path_Wvias,
          vdm = path_vdm, od = path_od, ibge = path_ibge)
.falta <- .req[!file.exists(.req)]
if (length(.falta)) {
  stop("Insumo(s) não encontrado(s):\n  ",
       paste(sprintf("[%s] %s", names(.falta), .falta), collapse = "\n  "),
       "\nAjuste DIR_DADOS / DIR_W_VIAS / DIR_CSVS no bloco de CONFIGURAÇÃO.")
}

# ==============================================================================
# PARTE 1: BASE - municípios, centroides, arestas de vizinhança (vias cruas)
# ==============================================================================
sp_municipios <- st_read(path_mun, quiet = TRUE) %>% mutate(id_inla = row_number())
n  <- nrow(sp_municipios)
cd <- as.character(sp_municipios$CD_MUN); nm <- as.character(sp_municipios$NM_MUN)
cd2id <- setNames(seq_len(n), cd)
coord_m <- st_coordinates(suppressWarnings(st_centroid(st_geometry(st_transform(sp_municipios, 31983)))))

mk <- function(a, b) paste(pmin(a, b), pmax(a, b), sep = "_")
nrm <- function(x){ x<-toupper(trimws(x)); x<-iconv(x,to="ASCII//TRANSLIT"); x<-gsub("[^A-Z ]","",x); trimws(gsub("\\s+"," ",x)) }

W_vias <- readRDS(path_Wvias); stopifnot(nrow(W_vias) == n)
Wt <- as(Matrix::triu(W_vias), "TsparseMatrix")
base_ed <- data.frame(iA = Wt@i + 1L, iB = Wt@j + 1L, vias = Wt@x)
base_ed <- base_ed[base_ed$iA != base_ed$iB & base_ed$vias > 0, ]
base_ed$key <- mk(base_ed$iA, base_ed$iB)
cat("Arestas base (vizinhança):", nrow(base_ed), "\n")

# ==============================================================================
# PARTE 2: CAMADAS QUE CRIAM ARESTAS NOVAS (proximidade e saltos OD)
# ==============================================================================
# --- proximidade kNN <= TETO_KM (não-base) ---
nb_k <- knn2nb(knearneigh(coord_m, k = K_VIZ), sym = TRUE)
prox <- do.call(rbind, lapply(seq_len(n), function(i){
  j <- nb_k[[i]]; j <- j[j > i]; if(!length(j)) return(NULL)
  d <- sqrt((coord_m[i,1]-coord_m[j,1])^2 + (coord_m[i,2]-coord_m[j,2])^2)/1000
  data.frame(iA=i, iB=j, dist_km=d)
}))
prox <- prox[prox$dist_km <= TETO_KM, ]; prox$key <- mk(prox$iA, prox$iB)
prox_nb <- prox[!(prox$key %in% base_ed$key), ]
prox_nb$proxval <- 1 / prox_nb$dist_km
cat("Atalhos de proximidade (não-base):", nrow(prox_nb), "\n")

# ==============================================================================
# PARTE 3: TABELAS DE FONTE POR PAR (cru), em id (iA<iB) -> key
# ==============================================================================
# --- OD (todos os pares; vira fonte e também cria saltos) ---
od_k <- data.frame(key=character(0), od=numeric(0))
tryCatch({
  mun_od <- st_read(path_od_mun, quiet = TRUE)
  dict <- setNames(as.character(mun_od$CD_IBGE), as.character(mun_od$NumeroMuni))
  od <- read.csv(path_od, stringsAsFactors = FALSE)
  od$idO <- cd2id[dict[as.character(od$MUNI_O)]]; od$idD <- cd2id[dict[as.character(od$MUNI_D)]]
  od <- od[!is.na(od$idO) & !is.na(od$idD) & od$idO != od$idD, ]
  od$key <- mk(od$idO, od$idD)
  od_k <- aggregate(viagens ~ key, od, sum); names(od_k)[2] <- "od"
}, error=function(e) cat("ERRO OD:", conditionMessage(e), "\n"))

# --- VDM (por nome -> cd -> id) ---
vdm_k <- data.frame(key=character(0), vdm=numeric(0))
tryCatch({
  vdm <- read.csv(path_vdm, stringsAsFactors = FALSE, fileEncoding = "UTF-8")
  vdm <- vdm[vdm$split_ok %in% c(TRUE,"TRUE") & !is.na(vdm$VDM2024), ]
  lut <- setNames(cd, nrm(nm))
  iA <- cd2id[lut[nrm(vdm$mun_A)]]; iB <- cd2id[lut[nrm(vdm$mun_B)]]
  ok <- !is.na(iA) & !is.na(iB) & iA != iB
  vdm_k <- aggregate(VDM2024 ~ key, data.frame(key=mk(iA[ok], iB[ok]), VDM2024=vdm$VDM2024[ok]), sum)
  names(vdm_k)[2] <- "vdm"
}, error=function(e) cat("ERRO VDM:", conditionMessage(e), "\n"))

# --- IBGE VAR07 (por código IBGE -> id) ---
ib_k <- data.frame(key=character(0), ibge=numeric(0))
tryCatch({
  ib <- read.csv(path_ibge, stringsAsFactors = FALSE, fileEncoding = "UTF-8")
  iA <- cd2id[as.character(ib$CODMUNDV_A)]; iB <- cd2id[as.character(ib$CODMUNDV_B)]
  ok <- !is.na(iA) & !is.na(iB) & iA != iB
  ib_k <- aggregate(VAR07 ~ key, data.frame(key=mk(iA[ok], iB[ok]), VAR07=ib$VAR07[ok]), sum)
  names(ib_k)[2] <- "ibge"
}, error=function(e) cat("ERRO IBGE:", conditionMessage(e), "\n"))

# ==============================================================================
# PARTE 4: UNIVERSO DE ARESTAS = base U proximidade U saltos-OD(não-base)
# ==============================================================================
od_nb <- od_k[!(od_k$key %in% base_ed$key), ]      # saltos OD: pares OD não-base
cat("Saltos de OD (não-base):", nrow(od_nb), "\n")

lookup <- rbind(
  base_ed[, c("key","iA","iB")],
  prox_nb[, c("key","iA","iB")],
  { p <- strsplit(od_nb$key, "_"); data.frame(key=od_nb$key,
                                              iA=as.integer(sapply(p,`[`,1)), iB=as.integer(sapply(p,`[`,2))) }
)
master <- lookup[!duplicated(lookup$key), ]
cat("TOTAL de arestas (Grafo #3):", nrow(master), "\n")

# ---- anexa valores CRUS de cada fonte (NA onde ausente) ----
master$vias <- base_ed$vias[match(master$key, base_ed$key)]            # só base
master$prox <- prox_nb$proxval[match(master$key, prox_nb$key)]         # só atalhos prox
master$od   <- od_k$od[match(master$key, od_k$key)]                    # qualquer aresta c/ OD
master$vdm  <- vdm_k$vdm[match(master$key, vdm_k$key)]                 # qualquer aresta c/ VDM
master$ibge <- ib_k$ibge[match(master$key, ib_k$key)]                  # qualquer aresta c/ IBGE

# ==============================================================================
# PARTE 5: RÉGUA ÚNICA - normaliza cada FONTE uma vez (global) e funde
# ==============================================================================
# z-score global por fonte: log1p + padronização sobre TODAS as arestas (não-NA)
zglobal <- function(x){
  ok <- !is.na(x); out <- rep(NA_real_, length(x))
  if (sum(ok) < 2) { out[ok] <- 0; return(out) }
  lx <- log1p(x[ok]); s <- sd(lx)
  out[ok] <- if (s > 0) (lx - mean(lx)) / s else 0
  out
}
Z <- cbind(vias = zglobal(master$vias), vdm = zglobal(master$vdm),
           od   = zglobal(master$od),   ibge = zglobal(master$ibge),
           prox = zglobal(master$prox))

# média ponderada das fontes PRESENTES (mesma régua p/ todas) -> 1 valor por aresta
pw <- pesos_fonte[colnames(Z)]
Wm <- matrix(pw, nrow(Z), ncol(Z), byrow = TRUE); Wm[is.na(Z)] <- NA
fundido <- rowSums(Z * Wm, na.rm = TRUE) / rowSums(Wm, na.rm = TRUE)
master$n_fontes <- rowSums(!is.na(Z))

# positivação ÚNICA p/ todas as arestas -> [floor_peso, 1]
rg <- range(fundido, na.rm = TRUE)
master$w <- if (diff(rg) > 0) floor_peso + (1-floor_peso)*(fundido-rg[1])/diff(rg) else 1

# rótulo de tipo (apenas p/ visualização)
is_base <- master$key %in% base_ed$key
is_prox <- master$key %in% prox_nb$key
is_odj  <- master$key %in% od_nb$key
master$tipo <- ifelse(is_base, "base",
                      ifelse(is_prox & is_odj, "proximidade+salto_od",
                             ifelse(is_prox, "proximidade", "salto_od")))

# ==============================================================================
# PARTE 6: MATRIZ W3, C = D - W3, restrições
# ==============================================================================
W3 <- sparseMatrix(i = master$iA, j = master$iB, x = master$w, dims = c(n,n)); W3 <- W3 + t(W3)
nb_chk <- mat2listw(as(W3 > 0, "dMatrix"), style="B", zero.policy=TRUE)$neighbours
comp <- n.comp.nb(nb_chk); ncomp <- comp$nc; compid <- comp$comp.id
stopifnot(sum(card(nb_chk) == 0) == 0)
D <- Diagonal(x = rowSums(W3)); C <- as(D - W3, "TsparseMatrix")
Aconstr <- matrix(0, nrow = ncomp, ncol = n)
for (cc in seq_len(ncomp)) Aconstr[cc, compid == cc] <- 1
econstr <- rep(0, ncomp)

# ==============================================================================
# PARTE 7: EXPORTAÇÃO  -> Resultados_Grafo3/
# ==============================================================================
saveRDS(W3,      file.path(DIR_RESULT, "W_grafo3.rds"))
saveRDS(C,       file.path(DIR_RESULT, "C_grafo3.rds"))
saveRDS(Aconstr, file.path(DIR_RESULT, "Aconstr_grafo3.rds"))
saveRDS(econstr, file.path(DIR_RESULT, "econstr_grafo3.rds"))
nb2INLA(file.path(DIR_RESULT, "sp_grafo3.graph"), nb_chk)
write.csv(data.frame(munA=nm[master$iA], munB=nm[master$iB],
                     cdA=cd[master$iA], cdB=cd[master$iB],
                     tipo=master$tipo, n_fontes=master$n_fontes,
                     vias=master$vias, vdm=master$vdm, od=master$od,
                     ibge=master$ibge, prox=master$prox, peso=master$w),
          file.path(DIR_RESULT, "grafo3_arestas.csv"), row.names = FALSE, fileEncoding = "UTF-8")

# ==============================================================================
# PARTE 8: VISUALIZAÇÃO  -> Resultados_Grafo3/
# ==============================================================================
coord_g <- st_coordinates(suppressWarnings(st_centroid(st_geometry(sp_municipios))))
edge_sf <- st_sf(tipo=factor(master$tipo), peso=master$w,
                 geometry = st_sfc(lapply(seq_len(nrow(master)), function(r)
                   st_linestring(rbind(coord_g[master$iA[r],], coord_g[master$iB[r],]))), crs = st_crs(sp_municipios)))
ordem <- c("base","proximidade","proximidade+salto_od","salto_od")
edge_sf$tipo <- factor(edge_sf$tipo, levels = intersect(ordem, levels(edge_sf$tipo)))
edge_sf <- edge_sf[order(edge_sf$tipo), ]
cores <- c("base"="grey78","proximidade"="darkorange","proximidade+salto_od"="purple","salto_od"="firebrick")

mapa_tipo <- ggplot() +
  geom_sf(data=sp_municipios, fill="gray98", color="gray88", linewidth=0.1) +
  geom_sf(data=edge_sf, aes(color=tipo), linewidth=0.45, alpha=0.8) +
  scale_color_manual(values=cores, name="tipo de aresta") +
  theme_minimal() +
  labs(title="Grafo #3 (régua única) - adjacência + atalhos",
       subtitle="Cinza=base | Laranja=proximidade<=150km | Vermelho=salto OD")
ggsave(file.path(DIR_RESULT,"Grafo3_mapa_tipos.png"), mapa_tipo, width=12, height=9, dpi=300, bg="white")

mapa_peso <- ggplot() +
  geom_sf(data=sp_municipios, fill="gray98", color="gray88", linewidth=0.1) +
  geom_sf(data=edge_sf[order(edge_sf$peso),], aes(color=peso, linewidth=peso), alpha=0.85) +
  scale_color_viridis_c(option="C", name="peso") + scale_linewidth(range=c(0.15,1.6), guide="none") +
  theme_minimal() +
  labs(title="Grafo #3 (régua única) - peso das arestas",
       subtitle="Pesos comparáveis: todas as fontes normalizadas na mesma régua")
ggsave(file.path(DIR_RESULT,"Grafo3_mapa_peso.png"), mapa_peso, width=12, height=9, dpi=300, bg="white")
print(mapa_tipo); print(mapa_peso)

# ==============================================================================
# PARTE 9: LOG  -> Resultados_Grafo3/
# ==============================================================================
sink(file.path(DIR_RESULT, "LOG_Grafo3.txt"))
cat("GRAFO #3 (RÉGUA ÚNICA) - LOG DE EXECUÇÃO\n")
cat("Gerado em:", as.character(Sys.time()), "\n\n")
cat("Config: TETO_KM=", TETO_KM, "| K_VIZ=", K_VIZ, "| floor_peso=", floor_peso, "\n")
cat("Pesos por fonte:", paste(names(pesos_fonte), pesos_fonte, sep="=", collapse=" "), "\n\n")
cat("NORMALIZAÇÃO: régua ÚNICA - cada fonte (vias, vdm, od, ibge, prox) é\n")
cat("  log1p+z-score sobre TODAS as arestas onde existe; fusão = média ponderada\n")
cat("  das fontes presentes; positivação única -> [", floor_peso, ", 1].\n\n")
cat("Arestas base (vizinhança)     :", nrow(base_ed), "\n")
cat("Atalhos de proximidade (novos):", nrow(prox_nb), "\n")
cat("Saltos de OD (novos)          :", nrow(od_nb), "\n")
cat("TOTAL de arestas              :", nrow(master), "\n\n")
cat("Distribuição por tipo:\n"); print(table(master$tipo))
cat("\nCobertura por fonte (nº de arestas com a fonte):\n")
for (f in c("vias","vdm","od","ibge","prox"))
  cat(sprintf("  %-5s: %d\n", f, sum(!is.na(master[[f]]))))
cat("\nDistribuição do nº de fontes por aresta:\n"); print(table(master$n_fontes))
cat("\nGrau por município:\n"); print(summary(card(nb_chk)))
cat("\nResumo do peso (régua única):\n"); print(summary(master$w))
cat("\nComponentes conexas:", ncomp, "\n")
sink()
cat("\nGrafo #3 (régua única) concluído. Arestas:", nrow(master), "| componentes:", ncomp, "\n")

# ==============================================================================
# PARTE 10: USO NO INLA (Caminho A) -- descomente
# ==============================================================================
# library(INLA)
# C_scaled <- inla.scale.model(C, constr = list(A = Aconstr, e = econstr))
# formula <- casos ~ 1 + f(id_inla, model="generic0", Cmatrix=C_scaled, constr=FALSE,
#   extraconstr=list(A=Aconstr, e=econstr), rankdef=ncomp,
#   hyper=list(prec=list(prior="pc.prec", param=c(1,0.01))))
# fit <- inla(formula, family="nbinomial", data=dados, E=dados$E,
#             control.compute=list(dic=TRUE, waic=TRUE)); summary(fit)