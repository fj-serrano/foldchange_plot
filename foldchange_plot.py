# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd 
import os
import sys
import argparse
from adjustText import adjust_text
import seaborn as sns

# FUNCIÓN QUE PERMITE AL USUARIO SELECCIONAR ARCHIVOS Y GUARDAR LAS RUTAS EN UNA LISTA
def select_archives():
    dic_files = {}
    estudios = []
    # Argumentos de entrada para el script que utiliza para localizar los archivos de entrada, 
    # el nombre del estudio, método de asignación de lecturas y localización de salida
    parser = argparse.ArgumentParser(description = "Procesa un archivo de entrada")
    parser.add_argument("-i", "--input", required = True, help = "Ruta del archivo de entrada")
    parser.add_argument("-o", "--output", required = True, help = "Ruta del archivo de salida")
    parser.add_argument("-s", "--study", help = "Elige el nombre del estudio", default = "SRP")
    parser.add_argument("-a", "--assignment", choices = ["de_rcadj","de_rcsa"], help = "Elige el metodo de asignacion", default = "de_rcadj")
    args = parser.parse_args()
    file_path = args.input
    output_path = args.output
    study = args.study
    assignment = args.assignment
    # Normaliza el path para evitar problemas con rutas relativas o absolutas
    file_path_clean = os.path.normpath(os.path.abspath(file_path))
    output_path = os.path.normpath(os.path.abspath(output_path))
    # Busca todas las carpetas cuyo nombre empiece por "study" y, dentro de cada una,
    # busca recursivamente una subcarpeta llamada igual que "assignment"
    # Si la encuentra, añade todas las .tsv que empiecen por DESeq2, edgeR o limma.
    metodos_prefijos = ("DESeq2", "edgeR", "limma")
    contador_carpetas = 1
    for root, dirs, files_in_dir in os.walk(file_path_clean):
        # Identifica directorios de estudio seleccionado
        if os.path.basename(root).startswith(study):
            # Si el estudio no está ya en la lista de estudios, añádelo
            if os.path.basename(root) not in estudios:
                estudios.append(os.path.basename(root))
            # Recorre las subcarpetas del estudio en busca de la carpeta de asignación
            found = False
            for subroot, subdirs, subfiles in os.walk(root):
                if os.path.basename(subroot) == assignment:
                    # Añade solo los archivos .tsv cuyos nombres empiecen por alguno de los metodos y en los que se enfrenten 2 condiciones (vs)
                    files = []
                    for file in subfiles:
                        if file.endswith(".tsv") and file.startswith(metodos_prefijos) and ("vs" in file):
                            file_path_full = os.path.join(subroot, file)
                            files.append(file_path_full)
                    dic_files[os.path.basename(root)] = files
                    # Al encontrar la carpeta de regresión dentro del estudio, no es necesario seguir buscando en ese estudio, se puede pasar al siguiente
                    print(f"{contador_carpetas}. Carpeta procesada: {os.path.basename(root)}")
                    contador_carpetas += 1
                    found = True
                    break
            # Si no se encontró la carpeta de regresión dentro del estudio, omitir
            if not found:
                continue
    return dic_files, output_path, assignment, estudios

# FUNCIÓN PARA DIAGRAMA DE BARRAS PARA EL PORCENTAJE DE miRNAs ÚNICOS EN UN INTERVALO DE FOLDCHANGE
def diagrama_barras(df_filt, name):
    sns.set_theme(style = "ticks")
    plt.rcParams.update({"font.size": 10})

    # Para el diagrama se definen los intervalos (bins) manteniendo simetria
    max_abs = np.ceil(np.max(np.abs(df_filt["log2FoldChange"])))
    # Para evitar posibles errores por falta de datos significativos, se establece un valor mínimo para max_abs
    if max_abs == 0:
        max_abs = 1
    intervalos = np.arange(-max_abs, max_abs + 0.5, 0.5)

    # Se asigna el tamaño del lienzo y se pinta el histograma con seaborn
    plt.figure(figsize = (7, 5))
    sns.histplot(data = df_filt, x = "log2FoldChange", bins = intervalos, stat = "probability", color = "#6EADFF", edgecolor = "white", alpha = 0.9)

    # Se representa una linea de simetría en x = 0 para facilitar la visualización de los cambios positivos y negativos
    plt.axvline(x = 0, color = "black", linestyle = "--", linewidth = 1.2, alpha = 0.7)
    plt.xlabel("Intervalo log2FoldChange")
    plt.ylabel("Porcentaje de miRNAs significativos")
    plt.title("Distribucion log2FoldChange para los miRNAs significativos", fontweight = "bold", pad = 12)
    sns.despine()
    plt.tight_layout()

    # Para guardar el archivo se le pedirá al usuario un nombre del estudio
    plt.savefig(f"log2FoldChange_0.05_{name}.png", dpi = 300, bbox_inches = "tight")
    plt.close()
    return

# FUNCIÓN PARA REPRESENTAR CONCENTRACIONES MEDIAS DE AMBAS CONDICIONES ENFRENTADAS
def distribucion_concentraciones(df, df_filt, name):
    plt.figure(figsize = (6.5, 6))
    columns = list(df.columns)

    # Se transforman las concentraciones a log10 para mejorar la visualización, pero solo se toman los valores mayores a 0
    x_all = np.log10(df[columns[1]].where(df[columns[1]] > 0))
    y_all = np.log10(df[columns[2]].where(df[columns[2]] > 0))
    x_sig = np.log10(df_filt[columns[1]].where(df_filt[columns[1]] > 0))
    y_sig = np.log10(df_filt[columns[2]].where(df_filt[columns[2]] > 0))

    # Se pintan los miRNAs no significativos en gris y los significativos en rojo
    plt.scatter(x_all, y_all, alpha = 0.3, color = "gray", s = 15, label = "No significativo")
    plt.scatter(x_sig, y_sig, alpha = 0.8, color = "orange", s = 15, label = "Significativo")

    # Se representan líneas de referencia en x = y para facilitar la visualización de los cambios entre las condiciones comparadas
    min_val = min(x_all.min(skipna = True), y_all.min(skipna = True))
    max_val = max(x_all.max(skipna = True), y_all.max(skipna = True))
    plt.plot([min_val, max_val], [min_val, max_val], color = "black", linestyle = "--", linewidth = 1.2, alpha = 0.7)

    # Se añaden etiquetas, título y leyenda, y se guarda la figura con el nombre del estudio proporcionado por el usuario
    label_x = columns[1].replace('_', ' ').title()
    label_y = columns[2].replace('_', ' ').title()

    plt.xlabel(f"log10({label_x})")
    plt.ylabel(f"log10({label_y})")
    plt.legend(loc = "upper left", fontsize = "9",  bbox_to_anchor = (1.02, 1), frameon = True)
    plt.title("Distribucion de Concentraciones", fontweight = "bold", pad = 12)
    sns.despine()
    plt.tight_layout()
    plt.savefig(f"distribucion_concentraciones_{name}.png", dpi = 300, bbox_inches = "tight")
    plt.close()
    return

# FUNCIÓN PARA REPRESENTAR GRÁFICO M-A, UN DIAGRAMA DE DISPERSIÓN PARA VER LA EXPRESIÓN DIFERENCIAL
def ma_plot(df, df_filt, name):
    plt.figure(figsize = (6.5, 6))
    df_plot = df.copy()
    df_filt_plot = df_filt.copy()
    columns = list(df_plot.columns)

    # Calculamos la media de las concentraciones para cada miRNA en ambas condiciones enfrentadas
    df_plot["mean_conc_mirnas"] = df[[columns[1], columns[2]]].mean(axis = 1)
    df_filt_plot["mean_conc_mirnas"] = df_filt[[columns[1], columns[2]]].mean(axis = 1)
    df_plot = df_plot.where(df_plot["mean_conc_mirnas"] > 0)
    df_filt_plot = df_filt_plot.where(df_filt_plot["mean_conc_mirnas"] > 0)

    # Y ahora se representa en el eje X la media de las concentraciones y en el eje Y el log2FoldChange
    x_all = np.log2(df_plot["mean_conc_mirnas"])
    y_all = df_plot["log2FoldChange"]
    x_sig_up = np.log2(df_filt_plot["mean_conc_mirnas"].where(df_filt_plot["log2FoldChange"] > 0))
    y_sig_up = df_filt_plot["log2FoldChange"].where(df_filt_plot["log2FoldChange"] > 0)
    x_sig_down = np.log2(df_filt_plot["mean_conc_mirnas"].where(df_filt_plot["log2FoldChange"] < 0))
    y_sig_down = df_filt_plot["log2FoldChange"].where(df_filt_plot["log2FoldChange"] < 0)

    # Se pintan los miRNAs no significativos en gris y los significativos en rojo
    plt.scatter(x_all, y_all, alpha = 0.3, color = "gray", s = 15, label = "No significativo")
    plt.scatter(x_sig_up, y_sig_up, alpha = 0.8, color = "orange", s = 15, label = "Sobreexpresado")
    plt.scatter(x_sig_down, y_sig_down, alpha = 0.8, color = "#6EADFF", s = 15, label = "Infraexpresado")

    # Se representan líneas de referencia en y = 0 para facilitar la visualización de los cambios entre las condiciones comparadas
    
    plt.axhline( y = 0, color = "black", linestyle = "--", linewidth = 1.2, alpha = 0.7)

    # Se añaden etiquetas, título y leyenda, y se guarda la figura con el nombre del estudio proporcionado por el usuario
    plt.xlabel(f"log2(Media Expresion)")
    plt.ylabel(f"log2FoldChange")
    plt.legend(loc = "upper left", fontsize = "9",  bbox_to_anchor = (1.02, 1), frameon = True)
    plt.title("MA plot", fontweight = "bold", pad = 12)
    sns.despine()
    plt.tight_layout()
    plt.savefig(f"ma_plot_{name}.png", dpi = 300, bbox_inches = "tight")
    plt.close()
    return

# FUNCIÓN PARA EL VOLCANO PLOT Y ARCHIVO CON INFORMACIÓN DE miRNAs SIGNIFICATIVOS
def volcano_plot(df, name):
    # Figura y ejes 
    fig, ax = plt.subplots(figsize = (12, 8))

    # Dibujo de los puntos no significativos en gris
    ax.scatter(df["log2FoldChange"], - np.log10(df["padj"]), alpha = 0.4, color = "gray", label = "No significativo")
    # Filtrado de significativos y dibujo de los puntos significativos en naranja y azul
    df_sig_sobre = df[(df["padj"] < 0.05) & (df["log2FoldChange"] > 1)].sort_values("padj")
    ax.scatter(df_sig_sobre["log2FoldChange"], - np.log10(df_sig_sobre["padj"]), alpha = 0.8, color = "orange", label = "Sobreexpresion")
    df_sig_infra = df[(df["padj"] < 0.05) & (df["log2FoldChange"] < -1)].sort_values("padj")
    ax.scatter(df_sig_infra["log2FoldChange"], - np.log10(df_sig_infra["padj"]), alpha = 0.8, color = "#6EADFF", label = "Infraexpresion")
    df_sig = df[(df["padj"] < 0.05) & (df["log2FoldChange"] > -1) & (df["log2FoldChange"] < 1)].sort_values("padj")
    ax.scatter(df_sig["log2FoldChange"], - np.log10(df_sig["padj"]), alpha = 0.8, color = "#31a354", label = "Estable")
    
    # Para añadir los nombres de los miRNAs con cambios significativos se filtra el data frame
    # para evitar una saturación de la imagen solo se mostrarán los nombres de los 5 que más han cambiado
    df_sig_total = df[(df["padj"] < 0.05) & (abs(df["log2FoldChange"]) > 1)].sort_values("padj")
    df_sig_volcano = df_sig_total.head(5)
    texts = []
    for index, row in df_sig_volcano.iterrows():
        x = row["log2FoldChange"]
        y = -np.log10(row["padj"])
        # pequeño desplazamiento inicial en y para evitar solapar exactamente sobre el punto
        txt = ax.text(x, y + 0.02, row["gene"], fontsize = 8, ha = "center", va = "bottom")
        texts.append(txt)

    # Ajustar las etiquetas usando la librería adjustText para evitar solapamientos
    # Se añade shrinkA/shrinkB y connectionstyle para evitar que las flechas atraviesen el texto
    adjust_text(
        texts,
        ax = ax,
        expand_text = (1.05, 1.2),
        expand_points = (1.2, 1.2),
        force_text = 0.6,
        force_points = 0.3,
        only_move = {"points": "xy", "texts": "xy"},
        # ShrinkA aumenta la separación del extremo en el texto (evita atravesarlo) y shrinkB la separación en el extremo del punto
        arrowprops = dict(arrowstyle = "->", color = "gray", lw = 0.5, shrinkA = 16, shrinkB = 2, connectionstyle  = "arc3,rad=0.15"),
        lim = 100
    )
 
    # Líneas de referencia
    ax.axvline(x = 1, color = "black")
    ax.axvline(x = -1, color = "black")
    ax.axhline(y = -np.log10(0.05), color = "black")

    # Etiquetas, título y leyenda
    ax.set_xlabel("log2FoldChange")
    ax.set_ylabel("-log10(padj)")
    ax.set_title("Volcano Plot", fontweight = "bold", pad = 12)
    ax.legend(loc = "upper left", fontsize = "9", bbox_to_anchor = (1.02, 1), frameon = True)
    fig.tight_layout()

    # Guardar el volcano plot con el nombre del estudio y mostrarlo
    fig.savefig(f"volcano_plot_{name}.png", dpi = 300, bbox_inches="tight")
    plt.close()
    return

def main():
    dic_files, output_path, assignment, estudios = select_archives()

    # Para cada archivo, se lee el data frame, se filtra por significancia y se representan los gráficos correspondientes
    for estudio in estudios:
        for file in dic_files[estudio]:
            df = pd.read_csv(file, sep = "\t")
            # Comprobamos antes de continuar que el data frame contiene las columnas log2FoldChange
            # y padj porque serán necesarias durante el análisis.
            if "log2FoldChange" not in df.columns or "padj" not in df.columns:
                print ("El archivo no es valido, debe contener log2FoldChange y padj.")
                sys.exit(2)

            name = os.path.basename(file).split(".")[0]
            # Cambiamos el directorio de trabajo a la ruta de salida proporcionada por el usuario
            output_path_new = output_path + f"/{estudio}/{assignment}"
            os.makedirs(output_path_new, exist_ok = True)
            os.chdir(output_path_new)

            # Eliminamos todos los posibles valores nulos
            df = df.dropna(subset = ["log2FoldChange", "padj"])

            # Ahora filtramos por significancia usando el valor padj que representa el valor p
            # ajustado y lo guardamos en la nueva variable df_filt. Además, solo tomamos aqullos tengan un log2FoldChange 
            # mayor a 1 o menor a -1 para centrarnos en los cambios más relevantes.
            df_filt = df[df["padj"] < 0.05]
            # Si no hay ningún miRNA significativo salta error
            if len(df_filt) == 0:
                print("Error: ningun miRNA significativo")
                continue
            else:
                # Guardar resultados de miRNAs significativos en un archivo .tsv con el nombre del estudio
                df_filt.to_csv(f"miRNAs_sig_{name}.tsv", sep = "\t", float_format = str, mode = "w", index = False)
                # Representación del diagrama de barras
                diagrama_barras(df_filt, name)
                # Representación de la distribución de concentraciones
                distribucion_concentraciones(df, df_filt, name)
                # Representación del volcano plot enfrentando FoldChange y pajus
                volcano_plot(df, name)
                # Representación del MA plot
                ma_plot(df, df_filt, name)


if __name__ == "__main__":
    main()