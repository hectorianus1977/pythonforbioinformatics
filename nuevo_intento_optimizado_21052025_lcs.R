# ============================
# CARGAR PAQUETES
# ============================
library(tidyverse)
library(Biostrings)
library(pheatmap)
library(ggplot2)
library(ggtext)
library(pbapply)

# ============================
# ARCHIVOS PERSONALIZADOS
# ============================
archivos_fasta <- c(
  "dsDNA.fasta",
  "ssRNA-RT.fasta",
  "ssRNA-.fasta",
  "ssRNA+.fasta",
  "ssDNA.fasta",
  "dsRNA.fasta",
  "dsDNA-RT.fasta"
)

# ============================
# FUNCIÓN PARA PROCESAR FASTA
# ============================
procesar_fasta <- function(archivo) {
  secuencias <- readAAStringSet(archivo)
  encabezados <- names(secuencias)
  secuencias_texto <- as.character(secuencias)
  
  tibble(Header = encabezados, Sequence = secuencias_texto) %>%
    mutate(
      ID = str_extract(Header, "^[^\\(\\s>]+"),
      Start = as.numeric(str_extract(Header, "(?<=\\()\\d+")),
      End = as.numeric(str_extract(Header, "(?<=-)\\d+(?=\\))")),
      MotifLength = End - Start + 1,
      ArchivoOrigen = basename(archivo),
      Categoria = str_replace(ArchivoOrigen, "-h\\.fasta$|\\.fasta$", "")
    ) %>%
    select(ID, Start, End, MotifLength, Sequence, ArchivoOrigen, Categoria)
}

# ============================
# PROCESAR TODOS LOS ARCHIVOS
# ============================
cat("📂 Cargando y procesando archivos FASTA...\n")
inicio <- Sys.time()

lista_resultados <- pblapply(archivos_fasta, procesar_fasta)
tabla_resultado <- bind_rows(lista_resultados)

cat("✅ Archivos FASTA procesados. Total secuencias: ", nrow(tabla_resultado), "\n")
cat("⏱️ Tiempo: ", round(Sys.time() - inicio, 2), " segundos\n\n")

write.csv(tabla_resultado, "secuencias_fusionadas_lcs.csv", row.names = FALSE)
cat("✅ Archivo fusionado guardado como: secuencias_fusionadas_lcs.csv\n")

# ============================
# FRECUENCIA DE AMINOÁCIDOS
# ============================
aminoacidos <- c("A","R","N","D","C","E","Q","G","H",
                 "I","L","K","M","F","P","S","T","W","Y","V")

cat("🧮 Calculando frecuencias de aminoácidos...\n")
inicio <- Sys.time()

calcular_frecuencias <- function(secuencia) {
  aa <- strsplit(secuencia, "")[[1]]
  freqs <- table(factor(aa, levels = aminoacidos))
  as.numeric(freqs)
}

frecuencia_aa <- pbapply::pblapply(tabla_resultado$Sequence, calcular_frecuencias)
frecuencia_aa <- do.call(rbind, frecuencia_aa)
colnames(frecuencia_aa) <- aminoacidos
tabla_frecuencias <- bind_cols(tabla_resultado, as.data.frame(frecuencia_aa))

cat("✅ Frecuencias calculadas.\n")
cat("⏱️ Tiempo: ", round(Sys.time() - inicio, 2), " segundos\n\n")

# ============================
# PORCENTAJE DE AMINOÁCIDOS POR SECUENCIA
# ============================
cat("📊 Calculando porcentaje de aminoácidos por secuencia...\n")
inicio <- Sys.time()

suma_aa <- rowSums(tabla_frecuencias[, aminoacidos])
tabla_porcentajes <- tabla_frecuencias
tabla_porcentajes[, aminoacidos] <- sweep(tabla_frecuencias[, aminoacidos], 1, suma_aa, FUN = "/") * 100

cat("✅ Porcentajes calculados.\n")
cat("⏱️ Tiempo: ", round(Sys.time() - inicio, 2), " segundos\n\n")

# ============================
# PROMEDIO DE PORCENTAJE POR CATEGORÍA
# ============================
cat("📈 Calculando promedios por categoría...\n")
inicio <- Sys.time()

promedios_categoria <- tabla_porcentajes %>%
  group_by(Categoria) %>%
  summarise(across(all_of(aminoacidos), \(x) mean(x, na.rm = TRUE))) %>%
  ungroup()

write.csv(promedios_categoria, "porcentajes_aminoacidos_por_categoria_lcs.csv", row.names = FALSE)
cat("✅ Promedios guardados en: porcentajes_aminoacidos_por_categoria_lcs.csv\n")
cat("⏱️ Tiempo: ", round(Sys.time() - inicio, 2), " segundos\n\n")

# ============================
# HEATMAP DE PORCENTAJE PROMEDIO
# ============================
cat("🧯 Generando heatmap de porcentajes promedio...\n")

matriz_heatmap <- as.data.frame(promedios_categoria)
rownames(matriz_heatmap) <- matriz_heatmap$Categoria
matriz_heatmap <- matriz_heatmap[, aminoacidos]

pheatmap(matriz_heatmap,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         display_numbers = TRUE,
         fontsize = 10,
         main = "🔥 Heatmap: Porcentaje promedio de aminoácidos por categoría")

png("heatmap_aminoacidos_por_categoria_lcs.png", width = 1000, height = 800)
pheatmap(matriz_heatmap,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         display_numbers = TRUE,
         fontsize = 10,
         main = "🔥 Heatmap: Porcentaje promedio de aminoácidos por categoría")


cat("✅ Heatmap guardado como: heatmap_aminoacidos_por_categoria_lcs.png\n\n")

# ============================
# ENRIQUECIMIENTO POR CATEGORÍA
# ============================
cat("🧪 Calculando enriquecimiento por categoría...\n")
inicio <- Sys.time()

promedio_global <- colMeans(matriz_heatmap)
enriquecimiento <- sweep(matriz_heatmap, 2, promedio_global, FUN = "-")

write.csv(enriquecimiento, "enriquecimiento_aminoacidos_por_categoria_lcs.csv")
cat("✅ Enriquecimiento guardado como: enriquecimiento_aminoacidos_por_categoria_lcs.csv\n")

png("enriquecimiento_heatmap_lcs.png", width = 1000, height = 800)
pheatmap(enriquecimiento,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         display_numbers = TRUE,
         fontsize = 10,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         main = "🎯 Enriquecimiento de aminoácidos por categoría (vs global)")

cat("✅ Heatmap de enriquecimiento guardado como: enriquecimiento_heatmap_lcs.png\n")
cat("⏱️ Tiempo: ", round(Sys.time() - inicio, 2), " segundos\n\n")

# ============================
# TOP 5 AMINOÁCIDOS POR CATEGORÍA
# ============================
cat("🏆 Calculando top 5 aminoácidos por categoría...\n")
inicio <- Sys.time()

top5_por_categoria <- promedios_categoria %>%
  pivot_longer(cols = all_of(aminoacidos), names_to = "Aminoacido", values_to = "Porcentaje") %>%
  group_by(Categoria) %>%
  arrange(desc(Porcentaje), .by_group = TRUE) %>%
  slice_head(n = 5) %>%
  ungroup()

write.csv(top5_por_categoria, "top5_aminoacidos_por_categoria_lcs.csv", row.names = FALSE)
cat("✅ Top 5 guardado como: top5_aminoacidos_por_categoria_lcs.csv\n")

ggplot(top5_por_categoria, aes(x = reorder(Aminoacido, Porcentaje), y = Porcentaje, fill = Aminoacido)) +
  geom_col(show.legend = FALSE) +
  geom_text(aes(label = paste0(round(Porcentaje, 1), "%")), hjust = -0.1, size = 3.5) +
  facet_wrap(~ Categoria, scales = "free_y") +
  coord_flip() +
  labs(
    title = "🎯 Top 5 aminoácidos más frecuentes por categoría_lcs",
    x = "Aminoácido",
    y = "Porcentaje (%)"
  ) +
  theme_minimal(base_size = 13) +
  theme(strip.text = element_text(face = "bold", size = 12))

ggsave("top5_aminoacidos_por_categoria_lcs.png", width = 10, height = 6)
cat("✅ Imagen guardada como: top5_aminoacidos_por_categoria_lcs.png\n")
