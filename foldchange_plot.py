# -*- coding: utf-8 -*-


import matplotlib.pyplot as plt
import numpy as np
# Pandas se usara para crear el data frame
import pandas as pd 
import sys
# Tkinter nos permite abrir el buscador de archivos para elegir el archivo a analizar y obtener la ruta
import tkinter 
from tkinter import filedialog

# Esta parte forma parte del import de tkinter 
tk = tkinter.Tk()
tk.withdraw()  

# FUNCIÓN PARA SELECCIONAR EL ARCHIVO CON LOS METADATOS
def selec_archive_to_df(): 
    # Guardamos en esta variable la ruta del excel seleccionado por el usuario
    file_path = filedialog.askopenfilename() 
    if not file_path:
        print("No se ha seleccionado ningun archivo.")
        sys.exit(1)
    # Creamos una variable vacía que contendrá la ruta sin / y con \
    file_path_clean = "" 

    for a in range(len(file_path)):
        if file_path[a] != "/":
            file_path_clean = file_path_clean + file_path[a]
        elif file_path[a] == "/":
            # Hace falta poner \\ porque sino lo interpreta como un carácter de escape
            file_path_clean = file_path_clean + "\\" 
    # Guardamos usando pandas la infomación del archivo .csv en un data frame      
    data_frame = pd.read_csv(file_path_clean, sep="\t") 
    return(data_frame)

# FUNCIÓN PARA DIAGRAMA DE BARRAS % DE miRNAs ÚNICOS VS INTERVALO DE FOLDCHANGE
def diagrama_barras(df_filt, name):
    # Generamos un intervalo (bins) en el que se clasificarán los datos a la hora de representarlos
    # con numpy hacemos el intervalo y con pd.cut dividimos los datos
    max_abs = np.ceil(np.max(np.abs(df_filt["log2FoldChange"])))
    interval = np.arange(-max_abs, max_abs  + 0.5, 0.5)
    categories = pd.cut(df_filt["log2FoldChange"], bins = interval)

    # Con value_counts de pandas obtenemos una serie de elementos de valores únicos
    # ordenados de forma que el primer elemento sea el que ocurre más frecuentemente
    # y excluyendo valores vacios. Además activamos la normalización para que los objetos 
    # contengan la frecuencia relativa. Y con sort_index ordenamos los objetos por label
    normal = categories.value_counts(normalize=True).sort_index()

    # Primero creamos una figura que mostrará el intervalo de log2FoldChange
    # contra el porcentaje de miRNAs significativos que muestran dicho valor
    plt.figure()
    normal.plot(kind="bar")

    plt.xlabel("Intervalo log2FoldChange")
    plt.ylabel("Proporcion de miRNAs significativos")
    plt.title("Distribucion de log2FoldChanbe para los miRNAs significativos (padj < 0.05)")
    plt.tight_layout()
    # Para guardar el archivo se le pedirá al usuario un nombre del estudio
    plt.savefig(f"log2FoldChanbe_0.05_{name}.png", bbox_inches="tight")
    plt.show()
    return

# FUNCIÓN PARA REPRESENTAR CONCENTRACIONES MEDIAS DE AMBAS CONDICIONES ENFRENTADAS
def distribucion_concentraciones(df, df_filt, name):
    plt.figure()
    columns = list(df.columns)
    plt.scatter(np.log10(df[columns[1]]), np.log10(df[columns[2]]), alpha = 0.4, color = "blue", label = "Azul (padj > 0.05)")
    plt.scatter(np.log10(df_filt[columns[1]]), np.log10(df_filt[columns[2]]), alpha = 0.4, color = "red", label = "Rojo (padj < 0.05)")
    plt.xlabel(f"log10({columns[1]})")
    plt.ylabel(f"log10({columns[2]})")
    plt.legend(loc="upper left")
    plt.title("Distribucion de Concentraciones")
    plt.tight_layout()
    plt.savefig(f"distribucion_concentraciones_{name}.png", bbox_inches="tight")
    plt.show()
    return

# FUNCIÓN PARA EL VOLCANO PLOT Y ARCHIVO CON INFORMACIÓN DE miRNAs SIGNIFICATIVOS
def volcano_plot(df, name):
    # Y luego creamos otra figura que mostrará los miRNAs de forma indivudal en un volcano plot
    # para ello representaremos el padj con el -log10 para aumentar la escala. Con alpha reducimos la opacidad de los puntos
    # facilitando la visualización
    plt.figure()
    plt.scatter(df["log2FoldChange"], -np.log10(df["padj"]), alpha = 0.4, color = "blue", label = "Azul (padj > 0.05)")

    # Para añadir los nombres de los miRNAs con cambios significativos se filtra el data frame
    # para evitar una saturación de la imagen solo se mostrarán los nombres de los 5 que más han cambiado
    df_sig_total = df[(df["padj"] < 0.05) & (abs(df["log2FoldChange"]) > 1)].sort_values("padj")
    df_sig_volcano = df[(df["padj"] < 0.05) & (abs(df["log2FoldChange"]) > 1)].sort_values("padj").head(5)
    plt.scatter(df_sig_total["log2FoldChange"], -np.log10(df_sig_total["padj"]), alpha = 0.4, color="red", label = "Rojo (padj < 0.05)")

    # Añado las etiquetas para los miRNAs significativos
    for index, row in df_sig_volcano.iterrows():
        plt.text(row["log2FoldChange"], -np.log10(row["padj"]), row["gene"], fontsize = 8)
    # Pinto lineas verticales en el gráfico que representan un cambio relevante biologicamente
    # en el que la expresión se duplica
    plt.axvline(x = 1, color = "black")
    plt.axvline(x =- 1, color = "black")
    # Pinto una linea horizontal en el gráfico que represnta el umbral de significancia
    plt.axhline(y = -np.log10(0.05), color = "black")
    plt.xlabel("log2FoldChange")
    plt.ylabel("-log10(padj)")
    plt.title("Volcano Plot")
    plt.tight_layout()
    plt.savefig(f"volcano_plot_{name}.png", bbox_inches="tight")
    plt.show()

    # Y guardo los valores significativos en un nuevo csv
    df_sig_total.to_csv(f"miRNAs_sig_{name}.tsv", sep="\t", float_format=str, mode="w", index=False)
    return

def main():
    df = selec_archive_to_df()

    # Comprobamos antes de continuar que el data frame contiene las columnas log2FoldChange
    # y padj porque serán necesarias durante el análisis.
    if "log2FoldChange" not in df.columns or "padj" not in df.columns:
        print ("El archivo no es valido, debe contener log2FoldChange y padj.")
        sys.exit(2)
    # Eliminamos todos los posibles valores nulos
    df = df.dropna(subset=["log2FoldChange", "padj"])

    # Ahora filtramos por significancia usando el valor padj que representa el valor p
    # ajustado y lo guardamos en la nueva variable df_filt
    df_filt = df[df["padj"] < 0.05]
    print (f"miRNAs unicos totales: {len(df)}")
    print (f"miRNAs unicos significativos (padj < 0.05): {len(df_filt)}")
    # Si no hay ningún miRNA significativo salta error, por tanto comprobar que haya
    if len(df_filt) == 0:
        print("Error: ningun miRNA significativo")
        sys.exit(3)
    # Usuario selecciona nombre de estudio
    name = input("Nombre del estudio: ")
    # Representación del diagrama de barras
    diagrama_barras(df_filt, name)
    # Representación de la distribución de concentraciones
    distribucion_concentraciones(df, df_filt, name)
    # Representación del volcano plot enfrentando FoldChange y pajus
    volcano_plot(df, name)
  

# Como buena práctica se hace esto antes de ejecutar main
if __name__ == "__main__":
    main()