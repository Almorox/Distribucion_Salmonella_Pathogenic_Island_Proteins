# Distribucion_Salmonella_Pathogenic_Island_Proteins

El script permite conocer la distribución de las proteínas que integran las trece Islas de Patogenicidad de Salmonella (Salmonella Pathogenic Islands, SPIs) en diferentes
proteomas pertenecientes a esta bacteria.

Toma como input dos archivos: 
- El archivo *query.fa*, con las proteínas cuya distribución se desea conocer.
- El archivo *proteomas_url*, con el nombre de las diferentes especies/subespecies de Salmonella 
y la url de NCBI para descargar sus respectivos proteomas.

Genera dos carpetas en el directorio del proyecto, una carpeta denominada data donde se almacenarán los archivos fasta descargados mediante wget, y los siguientes archivos:

- *proteome.fa* : contiene la concatenación de todos los archivos fasta. Es un multifasta para usarlo como subject en el análisis de blast.
- *Identifiers.tsv*: Es una copia del archivo ncbi_url_file
- *species_fasta_ID.tsv*: Es un archivo de dos columnas que contiene el nombre de cada especie junto a todos sus fasta_headers (IDs) de su respectivo proteome.fasta descargado.
- *query.fa*: es el archivo de query.fa utilizado como entrada, se copia a esta carpeta para que queden almacenadas y registradas las secuencias query utilizadas en el análisis.

Por otro lado los resultados se almacenan en la carpeta results que contiene, tras el análisis,
los siguientes archivos:

- *blast_result*: los resultados del blast en formato -outfmt"6 qseqid sseqid pident qcovs
sseq")
- *blast_result_filtered*: donde se almacenan solo aquellas líneas del blast_result que
pasen el cut off establecido de identidad y coverage
- *blast_result_final*: igual que blast_result_filtered, pero ahora incluyendo los nombres
de las especias a la que pertenece cada secuencia subject (hit) encontrada.

Por último, se crea un carpeta por cada proteína query, donde se almacenan los hits
encontrados para esa query en formato fasta (*proteina.fa*), su alineamiento con sus
proteínas homólogas (*proteina.aln*) y su árbol filogenético (*proteina.nw*). Estos
archivos tienen como identificador para cada secuencia no solo el header del fasta sino
también añadida la especie de Salmonella a la que pertenece.

También se crea un archivo log, donde se guardan las salidas del stdout y stderr, evitando
así el verbose de todos los programas usados mientras se ejecuta el script. Dicho archivo log se
almacena en la carpeta del proyecto (*output_folder*).

Además, el script tiene control sobre el número de argumentos (avisa al usuario si no ha introducido todos), avisa si el nombre del directorio declarado ya existe (le da al usuario la opción de abortar o de seguir con el proceso) y muestra finalmente en pantalla algunas estadísticas sobre los resultados obtenidos (número de hits y número de hits por query). Indica también dónde se han guardado los resultdos.
