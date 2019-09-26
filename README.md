#UNIX Assignment

##Data Inspection

###Attributes of `fang_et_al_genotypes`

```
head fang_et_al_genotypes.txt | awk -F "\t" '{print NF; exit}'
```

```
tail fang_et_al_genotypes.txt | awk -F "\t" '{print NF; exit}'
```

```
cut -f 1-10 fang_et_al_genotypes.txt | head | column -t
```

```
wc -l fang_et_al_genotypes.txt 
```

```
du -h fang_et_al_genotypes.txt
```

```
file fang_et_al_genotypes.txt
```

```
cut -f 4-986 fang_et_al_genotypes.txt | tail -n +2 | grep [^ATCG/] | less
```

From the inspection of this file, I learnt that:

* It has a total of 986 columns and 2783 rows
* Its size is approximately 6.1 millions of bytes
* It does not have any non ASCII characters
* It is a sequencing file of 2782 individuals (rows) for 983 snps (columns)
* It has some missing data represented by `?`


###Attributes of `snp_position.txt`

```
head snp_position.txt | awk -F "\t" '{print NF; exit}'
```

```
tail snp_position.txt | awk -F "\t" '{print NF; exit}'
```

```
wc -l snp_position.txt
```

```
head snp_position.txt | cut -f 1-10 | column -t
```

```
du -h snp_position.txt
```

```
file snp_position.txt 
```

By inspecting this file I learned that:

* It has 15 columns and 984 rows
* Its size is about 38,000 bytes
* It does not have any non ASCII characters
* It contains information about chromosome, gene and physical position of the 983 snps found in the previous file
* Unlike the previous file where each was a column, here each snp is a row

##Data Processing
###Processing of `snp_position.txt`

In order to merge the two files with the information needed, I used `cut`, `head`, `tail` and `sort` to format `snp_position.txt`. Below is the code used.
```
cut -f 1,3,4 snp_position.txt > snp_formatted_position.txt
```

```
(head -n 1 snp_formatted_position.txt && tail -n +2 snp_formatted_position.txt | sort -k1,1) > sorted_snp_formatted_position.txt
```

This code helped to select the columns snp_id, chromosome and position and sort the file by `snp_id` while keeping the header in the first row. 

I also did a post-processing inspection to check if the sorting was successful and any other errors. Below are the code lines.

```
tail -n +2 sorted_snp_formatted_position.txt  | sort -c -k1,1
```

```
head snp_formatted_position.txt | column -t
```


###Maize Data  

Coming the `fang_et_al_genotypes.txt`, I first had a look at how many distinct groups of individuals we have and their frequencies of appearance. 

```
sort -k3 fang_et_al_genotypes.txt | cut -f 3 | uniq -c | column -t
```

#####Filter, transpose and sort  

Next, I filtered out rows of the three target groups of maize, transposed and sorted the data by `snp_id` while always keeping the header in the first row.

```
awk -F "\t" '$3 ~ /ZMMIL|ZMMLR|ZMMMR|Group/' fang_et_al_genotypes.txt | cut -f 1,4-986 > maize_genotypes.txt
```

```
awk -f transpose.awk maize_genotypes.txt > transposed_maize_genotypes.txt
```

```
(head -n 1 transposed_maize_genotypes.txt && tail -n +2 transposed_maize_genotypes.txt | sort -k1,1) > sorted_transposed_maize_genotypes.txt
```

For quality check, I inspected both intermediate and final files as follows:

```
head maize_genotypes.txt | cut -f 1-10 | column -t
```

```
head maize_genotypes.txt | awk -F "\t" '{print NF; exit}'
```

```
cut -f 3 maize_genotypes.txt | sort | uniq -c
```

```
head transposed_maize_genotypes.txt | cut -f 1-10 | column -t 
```

```
head transposed_maize_genotypes.txt | awk -F "\t" '{print NF; exit}'
```

```
head sorted_transposed_maize_genotypes.txt | cut -f 1-10 | column -t 
```

```
tail -n +2 sorted_transposed_maize_genotypes.txt | sort -c -k1,1 
```
#####Merging snp and genotypes  

After all the files have been checked, I used `join` merged the snp_position and genotypes files for maize. Here, the code merged the files by their first columns while specifying that first row of each file is the header.

```
join --header -1 1 -2 1 -t $'\t' sorted_snp_formatted_position.txt sorted_transposed_maize_genotypes.txt > joined_maize_snp_and_genotypes.txt
```

#####Susbsetting per chromosome and position

Before subsetting the data per chromosme, I used `awk` to crop out all the rows with unknown and multiple position and keep only the snp with numeric position.

```
awk '$3 !~ /unknown|multiple/' joined_maize_snp_and_genotypes.txt > joined_maize_only_numeric_position.txt
```

I created a folder for maize and used `awk` and `for loop` to subset the data per chromosome both in increasing and decreasing orders as follows:

```
mkdir maize
```

```
for i in {1..10}; do (head -n 1 joined_maize_only_numeric_position.txt && awk '$2 == '$i'' joined_maize_only_numeric_position.txt | sort -k3,3n) > maize/maize_chr"$i"_increasing.txt; done
```

```
for i in {1..10}; do (head -n 1 joined_maize_only_numeric_position.txt && awk '$2 == '$i'' joined_maize_only_numeric_position.txt | sort -k 3nr | sed 's/?/-/g') > maize/maize_chr"$i"_decreasing.txt; done
```

```
awk '$3 ~ /unknown|Position/' joined_maize_snp_and_genotypes.txt > maize/maize_snp_with_unknown_position.txt
```

```
awk '$3 ~ /multiple|Position/' joined_maize_snp_and_genotypes.txt > maize/maize_snp_with_multiple_position.txt
```
#####Final inspection  

Finally, I did a quick inspection of all the 22 files generated for teosinte.

```
for i in {1..10}; do head -n 3 maize/maize_chr"$i"_increasing.txt | cut -f 1-7 | column -t; done
```

```
for i in {1..10}; do head -n 3 maize/maize_chr"$i"_decreasing.txt | cut -f 1-7 | column -t; done
```

```
head -n 3 maize/maize_snp_with_multiple_position.txt | cut -f 1-7 | column -t
```

```
head -n 3 maize/maize_snp_with_unknown_position.txt | cut -f 1-7 | column -t
```

###Teosinte Data

#####Filter, transposition and sort  

Just like for maize, I filtered out the three groups of interest in teosinte, transposed, sorted and inspected all the files.

```
awk -F "\t" '$3 ~ /ZMPBA|ZMPIL|ZMPJA|Group/' fang_et_al_genotypes.txt | cut -f 1,4-986 > teosinte_genotypes.txt
```

```
awk -f transpose.awk teosinte_genotypes.txt > transposed_teosinte_genotypes.txt
```

```
head transposed_teosinte_genotypes.txt | cut -f 1-10 | column -t
```

```
(head -n 1 transposed_teosinte_genotypes.txt && tail -n +2 transposed_teosinte_genotypes.txt | sort -k1,1) > sorted_transposed_teosinte_genotypes.txt
```

Below are the code lines for inspection

```
head teosinte_genotypes.txt | awk -F "\t" '{print NF; exit}
```

```
tail teosinte_genotypes.txt | awk -F "\t" '{print NF; exit}
```

```
cut -f 3 maize_genotypes.txt | sort | uniq -c
```

#####Merging the files  

Next, I merged and sorted the snp_position and genotypes files. 

```
join --header -1 1 -2 1 -t $'\t' sorted_snp_formatted_position.txt sorted_transposed_teosinte_genotypes.txt > joined_teosinte_snp_and_genotypes.txt
```

Similar to maize, I used `awk` to crop out all the rows with unknown and multiple position and remain only numeric position. Below is the code for that: 

```
awk '$3 !~ /unknown|multiple/' joined_teosinte_snp_and_genotypes.txt > joined_teosinte_only_numeric_position.txt
```

#####Subsetting files per chromosome and position  

Afterwards, I created a directory `teosinte` and used `awk` and `for loop` to subset the data per chromosome both in increasing and decreasing orders as follows:

```
mkdir teosinte
```

```
for i in {1..10}; do (head -n 1 joined_teosinte_only_numeric_position.txt && awk '$2 == '$i'' joined_teosinte_only_numeric_position.txt | sort -k3,3n) > teosinte/teosinte_chr"$i"_increasing.txt; done
```

```
for i in {1..10}; do (head -n 1 joined_teosinte_only_numeric_position.txt && awk '$2 == '$i'' joined_teosinte_only_numeric_position.txt | sort -k3,3nr | sed 's/?/-/g') > teosinte/teosinte_chr"$i"_decreasing.txt; done
```

Again, I used `awk` to subset the data for snp with unknown and multiple position.

```
awk '$3 ~ /unknown|Position/' joined_teosinte_snp_and_genotypes.txt > teosinte/teosinte_snp_with_unknown_position.txt
```

```
awk '$3 ~ /multiple|Position/' joined_teosinte_snp_and_genotypes.txt > teosinte/teosinte_snp_with_multiple_position.txt
```

#####Final inspection   

Finally, I did inspect all the 22 files generated for teosinte.

```
for i in {1..10}; do head -n 3 teosinte/teosinte_chr"$i"_increasing.txt | cut -f 1-7 | column -t; done
```

```
for i in {1..10}; do head -n 3 teosinte/teosinte_chr"$i"_decreasing.txt | cut -f 1-7 | column -t; done
```

```
head -n 3 teosinte/teosinte_snp_with_multiple_position.txt | cut -f 1-7 | column -t
```

```
head -n 3 teosinte/teosinte_snp_with_unknown_position.txt | cut -f 1-7 | column -t
```
