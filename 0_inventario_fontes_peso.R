# ==============================================================================
# INVENTÁRIO DAS FONTES DE PESO (para o GRAFO #2 ponderado por fluxo)
# Objetivo: LER e DESCREVER cada fonte (não monta o grafo ainda).
# Saídas (em Resultados_Inventario/):
#   - LOG_Inventario_Pesos.txt    (relatório de aprendizado)
#   - artesp_praca_diaria.csv     (tráfego médio diário por praça de pedágio)
#   - der_vdm_pares.csv           (pares de municípios + VDM 2024, do Excel DER)
#   - od_municipios_fluxo.csv     (viagens expandidas entre municípios, OD Metrô)
#   - ibge_ligacoes_sp.csv        (ligações de transporte entre municípios SP)
# ------------------------------------------------------------------------------
# CAMINHOS RELATIVOS: este script localiza a si mesmo e procura as pastas de
# dados AO LADO dele; toda a saída vai para "Resultados_Inventario/".
# ==============================================================================

library(sf)
library(dplyr)
library(readxl)
library(data.table)
library(stringr)
# foreign é usado p/ o .dbf da OD (carregado sob demanda na PARTE 3)

# ==============================================================================
# PARTE 0: ANCORAGEM DE CAMINHOS RELATIVOS
# ------------------------------------------------------------------------------
# Descobre a pasta onde ESTE script está (Rscript, source() ou botão Source do
# RStudio). As pastas de dados devem estar NESTA MESMA pasta. Se nada funcionar,
# usa o diretório de trabalho atual (getwd()).
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
DIR_DADOS  <- file.path(DIR_SCRIPT, "dados_entrada")   # pastas de dados ficam AQUI
DIR_RESULT <- file.path(DIR_SCRIPT, "Resultados_Inventario")
if (!dir.exists(DIR_RESULT)) dir.create(DIR_RESULT, recursive = TRUE)

# ------------------------------------------------------------------------------
# NOMES DAS PASTAS DE DADOS (dentro de dados_entrada/). Ajuste se renomear algo.
# ------------------------------------------------------------------------------
PASTA_ARTESP <- "ARTESP (Dados de Tráfego de Concessionárias  Pedágios)"
PASTA_VDM    <- "DER-SP (Volume Diário Médio - VDM)"
PASTA_OD     <- "Origem_Destinos_2023"
PASTA_IBGE   <- "Ligações Rodoviárias e Hidroviárias - Frequência de transporte público"

dir_artesp    <- file.path(DIR_DADOS, PASTA_ARTESP, "volume_trafego_diario_2024")
path_vdm_xlsx <- file.path(DIR_DADOS, PASTA_VDM, "VDM_2024.xlsx")
dir_vdm_shp   <- file.path(DIR_DADOS, PASTA_VDM, "Shapefiles - VDM 2024")
path_od_dbf   <- file.path(DIR_DADOS, PASTA_OD, "Site_190225_PesquisaOD2023", "Site_190225", "Banco2023_divulgacao_190225.dbf")
path_od_mun   <- file.path(DIR_DADOS, PASTA_OD, "Site_190225_PesquisaOD2023", "Site_190225", "002_Site Metro Mapas_190225", "Shape", "Municipios_2023.shp")
path_ibge     <- file.path(DIR_DADOS, PASTA_IBGE, "Base_de_dados_ligacoes_rodoviarias_e_hidroviarias_2016.xlsx")

path_log <- file.path(DIR_RESULT, "LOG_Inventario_Pesos.txt")

sf::sf_use_s2(FALSE)
linha <- function() cat(strrep("-", 72), "\n")

cat("Pasta do script :", DIR_SCRIPT, "\n")
cat("Pasta de saída  :", DIR_RESULT, "\n\n")

sink(path_log, split = TRUE)   # split=TRUE: imprime no console E no arquivo
cat("########################################################################\n")
cat("INVENTÁRIO DE FONTES DE PESO - GRAFO #2 (fluxo)\n")
cat("Gerado em:", as.character(Sys.time()), "\n")
cat("########################################################################\n\n")

# ==============================================================================
# PARTE 1: ARTESP - VOLUME DE PEDÁGIO DIÁRIO (20 lotes)
# Estrutura: 1 linha = praça x dia x sentido x tipo de pagamento.
# Colunas de veículo: tudo após as 7 primeiras (DATA..TP_PAGAMENTO).
# ATENÇÃO: o dado é por PRAÇA, não por par de municípios. Para virar peso de
#          aresta é preciso geolocalizar cada praça e dizer entre quais
#          municípios ela fica (passo manual/posterior).
# ==============================================================================
cat("[1] ARTESP - VOLUME DE PEDÁGIO\n"); linha()
tryCatch({
  arquivos <- list.files(dir_artesp, pattern = "VOLUME_PEDAGIADO.*\\.csv$", full.names = TRUE)
  cat("Arquivos encontrados:", length(arquivos), "\n")
  
  meta_um <- fread(arquivos[1], nrows = 5, encoding = "Latin-1")
  cat("Exemplo (", basename(arquivos[1]), ") - colunas:\n", sep = "")
  cat(paste(names(meta_um), collapse = ", "), "\n\n")
  
  fixas <- c("DATA","LOTE","PRACA","SENTIDO","TIPO_CABINE","TIPO_PASSAGEM","TP_PAGAMENTO")
  
  praca_tot <- rbindlist(lapply(arquivos, function(f) {
    dt <- fread(f, encoding = "Latin-1")
    veic <- setdiff(names(dt), fixas)
    veic <- veic[sapply(dt[, ..veic], is.numeric)]
    dt[, total_veic := rowSums(.SD, na.rm = TRUE), .SDcols = veic]
    # total no ano e nº de dias distintos -> média diária por praça
    dt[, .(veic_ano = sum(total_veic, na.rm = TRUE),
           n_dias   = uniqueN(DATA),
           lote     = LOTE[1]), by = PRACA]
  }))
  praca_tot[, veic_dia_medio := veic_ano / n_dias]
  
  cat("Total de praças distintas:", uniqueN(praca_tot$PRACA), "\n")
  cat("Tráfego médio diário por praça (top 10):\n")
  print(head(praca_tot[order(-veic_dia_medio), .(PRACA, lote, veic_dia_medio)], 10))
  fwrite(praca_tot, file.path(DIR_RESULT, "artesp_praca_diaria.csv"))
  cat("-> salvo: artesp_praca_diaria.csv\n")
}, error = function(e) cat("ERRO na leitura ARTESP:", conditionMessage(e), "\n"))
cat("\n")

# ==============================================================================
# PARTE 2: DER-SP - VDM (Volume Diário Médio)
# 2a) Excel: a coluna 'Trecho' já é "MunicipioA-MunicipioB" + VDM por ano.
# 2b) Shapefile: mesma info com GEOMETRIA (casável espacialmente).
# ==============================================================================
cat("[2] DER-SP - VDM\n"); linha()

# ---- 2a) Excel VDM ----
tryCatch({
  # cabeçalho real está na 2a linha (a 1a tem só os anos); leio sem nomes e atribuo por posição
  vdm <- read_excel(path_vdm_xlsx, sheet = "VDM", skip = 2, col_names = FALSE)
  nomes <- c("SP","Rodovia","Km","Trecho","Subtrecho","Km_ini","Km_fim",
             "Administrador","Data_concessao",
             "P2020","C2020","VDM2020","P2021","C2021","VDM2021",
             "P2022","C2022","VDM2022","P2023","C2023","VDM2023",
             "P2024","C2024","VDM2024")
  names(vdm)[seq_along(nomes)] <- nomes
  vdm <- vdm %>% filter(!is.na(Trecho))
  vdm$VDM2024 <- suppressWarnings(as.numeric(vdm$VDM2024))   # "-" vira NA
  
  # quebra 'Trecho' em par de municípios (CUIDADO: nomes com hífen, ex. Embu-Guaçu)
  partes <- str_split_fixed(vdm$Trecho, "-", 2)
  vdm$mun_A <- str_trim(partes[, 1])
  vdm$mun_B <- str_trim(partes[, 2])
  vdm$split_ok <- partes[, 2] != "" & !str_detect(partes[, 2], "-")
  
  cat("Linhas (subtrechos) no Excel VDM:", nrow(vdm), "\n")
  cat("Com VDM2024 preenchido:", sum(!is.na(vdm$VDM2024)), "\n")
  cat("Pares 'Trecho' que dividiram em 2 nomes limpos:", sum(vdm$split_ok), "\n")
  cat("Exemplos de Trecho->par + VDM2024:\n")
  print(head(vdm[, c("SP","Trecho","mun_A","mun_B","VDM2024")], 8))
  cat("\nAVISO: revise 'split_ok==FALSE' (nomes com hífen como Embu-Guaçu,\n")
  cat("       ou trechos com 3+ termos) antes de usar como chave de par.\n")
  write.csv(vdm[, c("SP","Rodovia","Trecho","mun_A","mun_B","VDM2024","split_ok")],
            file.path(DIR_RESULT, "der_vdm_pares.csv"), row.names = FALSE, fileEncoding = "UTF-8")
  cat("-> salvo: der_vdm_pares.csv\n")
}, error = function(e) cat("ERRO no Excel VDM:", conditionMessage(e), "\n"))
cat("\n")

# ---- 2b) Shapefiles VDM ----
tryCatch({
  shps <- list.files(dir_vdm_shp, pattern = "\\.shp$", full.names = TRUE)
  for (s in shps) {
    g <- st_read(s, quiet = TRUE)
    cat("Shapefile:", basename(s), "\n")
    cat("  feições:", nrow(g), "| geom:", as.character(unique(st_geometry_type(g))[1]),
        "| CRS:", st_crs(g)$input, "\n")
    cat("  colunas:", paste(setdiff(names(g), attr(g, "sf_column")), collapse = ", "), "\n")
    cat("  amostra de atributos:\n")
    print(head(st_drop_geometry(g), 3))
    cat("\n")
  }
}, error = function(e) cat("ERRO nos shapefiles VDM:", conditionMessage(e), "\n"))
cat("\n")

# ==============================================================================
# PARTE 3: OD METRÔ 2023 - FLUXO DE PESSOAS (resolve o Ferraz!)
# Banco...dbf = microdados de viagens. Agregamos viagens expandidas (FE_VIA)
# por par (MUNI_O, MUNI_D). Os códigos MUNI_* são da OD (não IBGE) -> mapear
# pelo shapefile Municipios_2023.
# ==============================================================================
cat("[3] OD METRÔ 2023 - FLUXO DE PESSOAS\n"); linha()
tryCatch({
  if (!requireNamespace("foreign", quietly = TRUE))
    stop("instale o pacote 'foreign' (install.packages('foreign'))")
  cat("Lendo .dbf da OD (grande, ~79MB; pode demorar)...\n")
  od <- foreign::read.dbf(path_od_dbf, as.is = TRUE)
  cat("Registros (viagens):", nrow(od), "| colunas:", ncol(od), "\n")
  cat("Colunas presentes (primeiras 60):\n")
  cat(paste(head(names(od), 60), collapse = ", "), "\n\n")
  
  # checa as colunas-chave esperadas
  chaves <- c("MUNI_O","MUNI_D","FE_VIA","MODOPRIN")
  cat("Colunas-chave presentes? ",
      paste(chaves, chaves %in% names(od), sep = "=", collapse = " | "), "\n")
  
  if (all(c("MUNI_O","MUNI_D","FE_VIA") %in% names(od))) {
    od$FE_VIA <- as.numeric(od$FE_VIA)
    fluxo <- od %>%
      filter(!is.na(MUNI_O), !is.na(MUNI_D), MUNI_O != MUNI_D) %>%
      group_by(MUNI_O, MUNI_D) %>%
      summarise(viagens = sum(FE_VIA, na.rm = TRUE), .groups = "drop") %>%
      arrange(desc(viagens))
    cat("Códigos de município de ORIGEM distintos:", length(unique(od$MUNI_O)), "\n")
    cat("Pares O-D distintos:", nrow(fluxo), "\n")
    cat("Top 10 pares por viagens expandidas/dia:\n"); print(head(fluxo, 10))
    write.csv(fluxo, file.path(DIR_RESULT, "od_municipios_fluxo.csv"), row.names = FALSE)
    cat("-> salvo: od_municipios_fluxo.csv\n")
  }
  rm(od); gc()
  
  # dicionário de códigos OD -> nome de município
  if (file.exists(path_od_mun)) {
    mun_od <- st_read(path_od_mun, quiet = TRUE)
    cat("\nMunicipios_2023 (dicionário OD) - colunas:",
        paste(setdiff(names(mun_od), attr(mun_od, "sf_column")), collapse = ", "), "\n")
    print(head(st_drop_geometry(mun_od), 5))
  }
}, error = function(e) cat("ERRO na OD:", conditionMessage(e), "\n"))
cat("\n")

# ==============================================================================
# PARTE 4: IBGE LIGAÇÕES RODOVIÁRIAS/HIDROVIÁRIAS 2016
# Leitura robusta: diagnostica o arquivo e tenta vários engines em cascata.
# ==============================================================================
cat("[4] IBGE - LIGAÇÕES RODOVIÁRIAS/HIDROVIÁRIAS 2016\n"); linha()
tryCatch({
  cat("=== DIAGNÓSTICO DO ARQUIVO ===\n")
  cat("Existe?      ", file.exists(path_ibge), "\n")
  if (!file.exists(path_ibge)) stop("Caminho não encontrado. Confira o nome da pasta/arquivo.")
  cat("Tamanho (MB):", round(file.size(path_ibge) / 1e6, 2),
      " (esperado ~9.2 MB; muito menor = download truncado)\n")
  
  # Um .xlsx é um ZIP. Os 2 primeiros bytes devem ser "PK" (0x50 0x4B).
  mb <- readBin(path_ibge, what = "raw", n = 4)
  assinatura <- rawToChar(mb[1:2])
  cat("Assinatura inicial:", assinatura,
      if (assinatura == "PK") " (OK, é um zip/xlsx válido)\n"
      else " (NÃO é 'PK' -> não é xlsx válido: pode ser .xls antigo, .csv renomeado, ou corrompido)\n")
  cat("\n")
  
  # Leitura em cascata: para no primeiro método que funcionar.
  ler_ibge <- function(path) {
    tentativas <- list(
      "readxl::read_excel (auto)" = function() readxl::read_excel(path, guess_max = 100000),
      "readxl + format xlsx"      = function() readxl::read_xlsx(path, guess_max = 100000),
      "openxlsx::read.xlsx"       = function() {
        if (!requireNamespace("openxlsx", quietly = TRUE)) stop("pacote openxlsx ausente")
        openxlsx::read.xlsx(path) },
      "readxl como .xls (xls disfarçado)" = function() readxl::read_xls(path),
      "cópia local + readxl"      = function() {
        tmp <- file.path(tempdir(), "ibge_lig.xlsx")
        file.copy(path, tmp, overwrite = TRUE)   # contorna lock do OneDrive/Excel
        readxl::read_excel(tmp, guess_max = 100000) }
    )
    for (nm in names(tentativas)) {
      res <- try(tentativas[[nm]](), silent = TRUE)
      if (!inherits(res, "try-error") && !is.null(res) && nrow(res) > 0) {
        cat(">> SUCESSO via:", nm, "\n\n"); return(res)
      } else {
        cat("   falhou:", nm,
            if (inherits(res, "try-error")) paste0(" [", conditionMessage(attr(res, "condition")), "]"),
            "\n")
      }
    }
    stop("Nenhum método conseguiu ler o arquivo. Veja o diagnóstico acima.")
  }
  
  lig <- ler_ibge(path_ibge)
  
  cat("Linhas:", nrow(lig), "| colunas:", ncol(lig), "\n")
  cat("Nomes das colunas:\n"); cat(paste(names(lig), collapse = ", "), "\n\n")
  
  vars <- grep("^VAR", names(lig), value = TRUE)
  cat("Tipos das colunas VAR* (p/ achar frequência/tempo/custo):\n")
  for (v in vars) cat("  ", v, ":", class(lig[[v]])[1],
                      "| ex:", paste(head(unique(lig[[v]]), 3), collapse = " / "), "\n")
  
  cat("\nAmostra (pares + primeiras VARs):\n")
  cols_show <- intersect(c("UF_A","NOMEMUN_A","UF_B","NOMEMUN_B", head(vars, 6)), names(lig))
  print(head(lig[, cols_show], 6))
  
  # Recorte SP-SP
  if (all(c("UF_A","UF_B") %in% names(lig))) {
    lig_sp <- lig %>% filter(UF_A == "SP", UF_B == "SP")
    cat("\nLigações SP-SP:", nrow(lig_sp), "\n")
    write.csv(lig_sp, file.path(DIR_RESULT, "ibge_ligacoes_sp.csv"),
              row.names = FALSE, fileEncoding = "UTF-8")
    cat("-> salvo: ibge_ligacoes_sp.csv\n")
  } else {
    cat("\nAVISO: colunas UF_A/UF_B não encontradas com esse nome; ajuste o filtro SP.\n")
  }
  
  cat("\nNOTA (dicionário IBGE confirmado): VAR03=custo(R$), VAR04=tempo(min),\n")
  cat("     VAR05=freq.hidroviária, VAR06=freq.rodoviária, VAR07=FREQUÊNCIA TOTAL\n")
  cat("     de saídas (esta é a usada como peso), VAR08-11=coordenadas, VAR14=custo/tempo.\n")
}, error = function(e) cat("ERRO no Excel IBGE:", conditionMessage(e), "\n"))
cat("\n")

# ==============================================================================
# PARTE 5: SÍNTESE / PRÓXIMO PASSO
# ==============================================================================
cat("[5] SÍNTESE - COMO CADA FONTE VIRA PESO DE ARESTA (W) NO GRAFO #2\n"); linha()
cat("FONTE              CHAVE DE LIGAÇÃO              ESFORÇO   COBERTURA\n")
cat("DER VDM (xlsx/shp) coluna 'Trecho' = par A-B    BAIXO     rodovias DER (estado)\n")
cat("IBGE Ligações     CODMUNDV_A / CODMUNDV_B        BAIXO     todo o estado (2016)\n")
cat("OD Metrô          MUNI_O / MUNI_D (cód OD)       MÉDIO     RMSP (resolve Ferraz)\n")
cat("ARTESP            praça -> precisa geolocalizar  ALTO      eixos concessionados\n")
cat("\nMais fáceis de plugar primeiro: DER VDM e IBGE Ligações (já são por par).\n")
cat("OD entra para densificar a RMSP. ARTESP por último (exige geocodificar praças).\n")
sink()

cat("\nInventário concluído. Resultados em:\n", DIR_RESULT, "\n")