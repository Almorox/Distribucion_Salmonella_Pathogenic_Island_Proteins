#!/bin/bash
#Ejercicio Evaluable bash v1.0
#Lucía Almorox Antón 2021

# This script allow to performs blastp analysis using a urls_file to download the subject proteomes
# this file have to contain two columns, first column contains the species (identifier) for that proteome
# and the url. The script will generate a output project folder which contains two subfolders, data (here
# downloaded proteomes are stored and input fasta files) and results (here blast result, in addition a folder is created
# for every query protein to store blast result, aligment and trees for each protein)


###Función para imprimir ayuda en la terminal
ayuda () {
   echo 'Ejercicio Evaluable bash v1.0'
   echo 'Lucía Almorox Antón 2021'
   echo -e '\nusage: ejercicio.sh <query_sequences.fa> <ncbi_urls_file> <output_folder> <blast-identity> <blast-coverage>\n'
   echo -e 'query_sequences.fa : a fasta/multifasta file containing protein query sequences for blast'
   echo -e 'ncbi_urls_file     : a text plain file containing species name and url to download fasta protein file'
   echo -e 'output_folder          : folder in which data and results will be stored'
   echo -e 'blast-identity     : sequence identity cut off value 0-100'
   echo -e 'blast-coverage     : sequence coverage cut off value 0-100\n'
}

###Control del mensaje de ayuda
#Si el primer argumento es -h o -help, se invoca a la función ayuda y se corta la ejecución del script con número de terminación cero.
case $1 in
"-h") ayuda
      exit 0 ;;
"-help") ayuda
         exit 0 ;;
*);;
esac

###Control del número de argumentos con el que se invoca al script.
#Si el número de argumentos es menor de cinco, la terminal avisa del error, invocamos al función ayuda y terminamos la ejecución del script.
#número de terminación distinto de cero porque ha habido un error.
if [ $# -lt 5 ]
then 
echo "error:too few arguments"
ayuda
exit 1
fi

#Asignamos un nombre a cada argumento
query=$1
ncbi_urls_file=$2
project_name=$3
iden=$4 #70
cov=$5 #40

#Si el nombre del proyecto no existe, creamos una carpeta con su nombre y todos los directorios y ficheros donde guardaremos los outputs.
#Si sí que exite, la terminal avisará del error y se cortará la ejecución del script (número de terminación distinto de cero).
if test ! -d $project_name
then
  mkdir $project_name
  mkdir $project_name/Data
  mkdir $project_name/Results
  touch $project_name/log

  touch $project_name/Data/Identifiers.tsv
  cp $query $project_name/Data       #copiamos el archivo del primer parámetro a esta nueva dirección. 
  
  if [[ $query != query.fa ]]   #Si dicho archivo no se llama query.fa, le pondremos ese nombre
  then
    mv $project_name/Data/$query $project_name/Data/query.fa
  fi
 
  touch $project_name/Data/species_fasta_ID.tsv
  touch $project_name/Data/subject_proteome.fa

  touch $project_name/Results/blast_result.txt
  touch $project_name/Results/blast_result_filtered.tsv
  touch $project_name/Results/blast_result_final.tsv

  cat $project_name/Data/query.fa | grep ">" | sed 's/>//' | while read line  #Creamos una carpeta para cada proteína query
  do
    mkdir $project_name/Results/$line
  done

else
  print "Project name ($project_name) already exists as a directory. Script execution is going to be inturrepted."
  exit 1
fi

###Parseo de $ncbi_urls_file y creación del contenido de proteome.fa, species_fasta.tsv, Identifiers.tsv
#Los nombres de las especies quedan guardadas en el archivo temporal temp.$$ 
#guardo el contenido de $ncbi_urls_file separado por espacios para que sea más fácil el grep que hago abajo, en el que contemplo que haya un espacio 
#(que en el archivo original es un tab) para tener en cuenta dónde acaba el nombre de la especie (si no hiciera eso, al buscar la especie Tiphy, Tiphymurium_LT2
#también haría match).
echo -e "$0 is Running...\n"
echo "Fasta files will be stored in /$project_name/data/ "
echo "Generating species_fasta_ID.tsv file.. "

cut -f1 $ncbi_urls_file >> temp.$$
cat $ncbi_urls_file | tr "\t" " " > newurls.$$
cat temp.$$ | while read line #para cada especie ...
do 
  esp=$line
  link=$(grep "$esp " newurls.$$) #$link almacena la línea que contiene tanto el nombre de la especie como el link de descarga
  cd $project_name                #se descarga y descomprime en este directorio
  wget $link >> log 2>&1 
  gzip -d *.gz
  for file in $(ls *.faa)               
  do
     mv $file Data                #se mueve el archivo a Data, de forma que la siguiente vez que recorra el bucle, en "for file in $(ls)" ya no aparecerá este archivo (solo aparecerá el correspondiente a la especie con la que está el bucle en ese momento).
     cat Data/$file >> Data/subject_proteome.fa
     echo "$project_name/Data/$file" >> Data/temp2.$$     #el nombre del archivo (GCF...) queda almacenado en un fichero temporal (va precedido por $project_name/Data)(lo necesitaré para el archivo Identifiers.tsv)
     awk -v essp="$esp" '/>lcl/ {print essp, $1}' Data/$file >> Data/species_fasta_ID.tsv  #se almacenan todos los IDs del archivo precedidos por el nombre de la especie.
  done
  cd ..   #después de hacer esto para una especie, se vuelve al directorio en el que se está ejecutando el script
done 

rm newurls.$$
paste temp.$$ $project_name/Data/temp2.$$ > $project_name/Data/Identifiers.tsv  #en una columna quedan los nombres de especies y en la segunda $project_name/Data/GCF...
rm temp.$$
rm $project_name/Data/temp2.$$


#### BLAST SECTION ####
#Análisis con blast.
echo "Running Blast analsys and result filtering... "

blastp -query $project_name/Data/query.fa -subject $project_name/Data/subject_proteome.fa -outfmt "6 qseqid sseqid pident qcovs sseq" > $project_name/Results/blast_result.txt 2>>log

#Filtrado del resultado, en base a los valores de identity y coverage dados al ejecutar el script.
cd $project_name/Results
awk -v ident="$iden" -v cover="$cov" '$3 >= ident && $4 >= cover {print}' blast_result.txt  >> blast_result_filtered.tsv
cd ..
cd ..

#Creación de blast_result_final

#pim.$$ almacena la primera columna del archivo filtrado, pam.$$ el resto. pum$$ almacena la columna de headers: 
#recorremos este archivo para buscar cada header en el archivo $project_name/Data/species_fasta_ID.tsv 
#y con cada match nos llevamos solo el nombre de la especie, que queda almacenado en especitas.$$ 
#(tendremos al final, en ese archivo temporal, una columna con las especies correspondientes a cada header que 
#aparece en el archivo filtrado de blast). Mediante pastes quedan pim.$$, especitas$$ y pam$$ unidos. 

awk '{print $1}' $project_name/Results/blast_result_filtered.tsv >> pim.$$
awk '{print $2 "\t" $3 "\t" $4 "\t" $5}' $project_name/Results/blast_result_filtered.tsv >> pam.$$
awk '{print $2}' $project_name/Results/blast_result_filtered.tsv >> pum.$$
cat pum.$$ | while read line
do
  code=$line
  maa=$(grep -m 1 "$code" $project_name/Data/species_fasta_ID.tsv) #solo queremos el nombre de la especie una vez, es decir, solo el primer match, por eso -m 1
  echo $maa | cut -d ' ' -f1 >> especitas.$$
done

paste pim.$$ especitas.$$ >> casicasi.$$
paste casicasi.$$ pam.$$ >> $project_name/Results/blast_result_final.tsv

rm pim.$$
rm pam.$$
rm pum.$$
rm especitas.$$
rm casicasi.$$

#Creación del contenido para cada query
#almacenamos las proteínas query, de forma única, en una lista que recorremos para buscar cada una de ellas en blast_result_final.tsv, de donde, con cada match, 
#nos llevamos: nombre de la especie, header y  ecuencia. Todo ello se redirecciona a un nuevo archivo .fa que se crea en la carpeta para esa proteína query.
#me parece interesante comentar que no pude hacer /quer/ para buscar los matches, porque las barras hacen que "quer" se interprete literal
#así que $0~quer es una alternativa. 
awk '{print $1}' $project_name/Results/blast_result_final.tsv | uniq >> $project_name/Results/uniq_query_list.txt
cat $project_name/Results/uniq_query_list.txt | while read line
do
 queryy=$line
 awk -v quer="$queryy" '$0~quer{print ">" $2 "_" $3 "\n" $6}' $project_name/Results/blast_result_final.tsv  >> $project_name/Results/$queryy/$queryy.fa
done


### MUSCLE SECTION ####
#aprovecho esta sección en la que recorro los archivos .fa de cada query (con los que hago el alineamiento y el árbol) 
#para contar, además, los hits de cada proteína y los hits totales.
echo "Making MUSCLE alignment and phylogenetic trees..."
cd $project_name/Results
tot_hits=0
for dir in $(ls -d */)
do
  cd $dir          #entramos en el directorio de una proteína query
  for hitfafile in $(ls *.fa)
  do
    prot=$(echo $hitfafile | cut -d '.' -f1)  #si al fichero SiiA.fa, por ejemplo, le quitamos .fa, nos quedamos con la proteína SiiA
    num_hits=$(grep -c ">" $hitfafile)
    tot_hits=`expr $tot_hits + $num_hits`   #sumamos los hits de esa proteína al contador de hits totales
    echo "$num_hits $prot " >> ../../final_terminal.$$   #en este fichero temporal vamos almacenando cada nombre de proteína query precedido por su número de hits
    muscle -in $hitfafile -out $prot.aln >> ../../log 2>&1
    muscle -maketree -in $prot.aln -out $prot.aln.nw -cluster neighborjoining >> ../../log 2>&1   
  done
  cd ..  #volvemos al directorio $project_name/Results
done
cd ..
cd ..    #volvemos al directorio desde donde hemos ejecutado el script (donde se crea el proyecto)

### TERMINAL PRINTING ####

echo -e "*** Done!! ****\n"
echo "Results are available at /$project_name/Results/"
echo "Total Blast hits  $tot_hits"
echo "hits were found for query proteins:"
cat $project_name/final_terminal.$$ | tr "\n" " "
rm $project_name/final_terminal.$$
