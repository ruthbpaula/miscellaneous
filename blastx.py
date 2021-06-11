import os
import subprocess

a=["cbrongniartii_completo","cbrongniartii_completo","cbrongniartii_completo","cmilitaris_completo","cmilitaris_completo","cmilitaris_completo","efestucae_nuclear","efestucae_nuclear","efestucae_nuclear","fcircinatum_nuclear","fcircinatum_nuclear","fcircinatum_nuclear","fgraminearum_nuclear","fgraminearum_nuclear","fgraminearum_nuclear","fsolani_nuclear","fsolani_nuclear","fsolani_nuclear","fverticillioides_nuclear","fverticillioides_nuclear","fverticillioides_nuclear","manisopliae_completo","manisopliae_completo","manisopliae_completo","mrobertsii_completo","mrobertsii_completo","mrobertsii_completo","tasperellum_nuclear","tasperellum_nuclear","tasperellum_nuclear","tatroviride_nuclear","tatroviride_nuclear","tatroviride_nuclear","tgamsii_nuclear","tgamsii_nuclear","tgamsii_nuclear","thamatum_completo","thamatum_completo","thamatum_completo","treesei_nuclear","treesei_nuclear","treesei_nuclear"]
b=["LAGLIDADG","GIY","YIG","LAGLIDADG","GIY","YIG","LAGLIDADG","GIY","YIG","LAGLIDADG","GIY","YIG","LAGLIDADG","GIY","YIG","LAGLIDADG","GIY","YIG","LAGLIDADG","GIY","YIG","LAGLIDADG","GIY","YIG","LAGLIDADG","GIY","YIG","LAGLIDADG","GIY","YIG","LAGLIDADG","GIY","YIG","LAGLIDADG","GIY","YIG","LAGLIDADG","GIY","YIG","LAGLIDADG","GIY","YIG"]
c=["cbro","cbro","cbro","cmil","cmil","cmil","efes","efes","efes","fcir","fcir","fcir","fgra","fgra","fgra","fsol","fsol","fsol","fver","fver","fver","mani","mani","mani","mrob","mrob","mrob","tasp","tasp","tasp","tatr","tatr","tatr","tgam","tgam","tgam","tham","tham","tham","tree","tree","tree"]

for i in range(len(a)):
        cmd="/mnt/c/ncbi-blast-2.9.0+/bin/blastx.exe -query "+a[i]+".fasta -db "+b[i]+".fasta -num_threads 4 -outfmt 6 -out blast_"+c[i]+"_"+b[i]+".txt"
        print(cmd)
        process = subprocess.Popen(cmd, shell=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        # wait for the process to terminate
        out, err = process.communicate()
        errcode = process.returncode
        print(out)
        print(err)
