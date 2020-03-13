############################################
# shell script to run an entire simulation #
############################################

# 23-Feb-2020
# This shell script was written for 
# Neighbor GWAS: incorporating neighbor genotypic identity into genome-wide association studies of field herbivory on Arabidopsis thaliana
# which was co-authored by Yasuhiro Sato, Eiji Yamamoto, Kentaro K. Shimizu and Atsushi J. Nagano
# correspondence: sato.yasuhiro.36c@kyoto-u.jp

Rscript neighborGWASsimul_prepGeno.R

#run at 300 SNPs
Rscript neighborGWAS250kSimulation.R 3 300 6 3 1 0.6 0.8 1&
Rscript neighborGWAS250kSimulation.R 3 300 6 3 1 0.3 0.8 2&
Rscript neighborGWAS250kSimulation.R 3 300 6 3 1 0.3 0.4 3&
Rscript neighborGWAS250kSimulation.R 3 300 6 3 1 0.1 0.4 4&

Rscript neighborGWAS250kSimulation.R 1 300 6 3 1 0.6 0.8 5&
Rscript neighborGWAS250kSimulation.R 1 300 6 3 1 0.3 0.8 6&
Rscript neighborGWAS250kSimulation.R 1 300 6 3 1 0.3 0.4 7&
Rscript neighborGWAS250kSimulation.R 1 300 6 3 1 0.1 0.4 8&

Rscript neighborGWAS250kSimulation.R 0.25 300 6 3 1 0.6 0.8 9&
Rscript neighborGWAS250kSimulation.R 0.25 300 6 3 1 0.3 0.8 10&
Rscript neighborGWAS250kSimulation.R 0.25 300 6 3 1 0.3 0.4 11&
Rscript neighborGWAS250kSimulation.R 0.25 300 6 3 1 0.1 0.4 12&

Rscript neighborGWAS250kSimulation.R 3 300 4 4 1 0.6 0.8 13&
Rscript neighborGWAS250kSimulation.R 3 300 4 4 1 0.3 0.8 14&
Rscript neighborGWAS250kSimulation.R 3 300 4 4 1 0.3 0.4 15&
Rscript neighborGWAS250kSimulation.R 3 300 4 4 1 0.1 0.4 16&

Rscript neighborGWAS250kSimulation.R 1 300 4 4 1 0.6 0.8 17&
Rscript neighborGWAS250kSimulation.R 1 300 4 4 1 0.3 0.8 18&
Rscript neighborGWAS250kSimulation.R 1 300 4 4 1 0.3 0.4 19&
Rscript neighborGWAS250kSimulation.R 1 300 4 4 1 0.1 0.4 20&

Rscript neighborGWAS250kSimulation.R 0.25 300 4 4 1 0.6 0.8 21&
Rscript neighborGWAS250kSimulation.R 0.25 300 4 4 1 0.3 0.8 22&
Rscript neighborGWAS250kSimulation.R 0.25 300 4 4 1 0.3 0.4 23&
Rscript neighborGWAS250kSimulation.R 0.25 300 4 4 1 0.1 0.4 24&

Rscript neighborGWAS250kSimulation.R 3 300 3 6 1 0.6 0.8 25&
Rscript neighborGWAS250kSimulation.R 3 300 3 6 1 0.3 0.8 26&
Rscript neighborGWAS250kSimulation.R 3 300 3 6 1 0.3 0.4 27&
wait
Rscript neighborGWAS250kSimulation.R 3 300 3 6 1 0.1 0.4 28&

Rscript neighborGWAS250kSimulation.R 1 300 3 6 1 0.6 0.8 29&
Rscript neighborGWAS250kSimulation.R 1 300 3 6 1 0.3 0.8 30&
Rscript neighborGWAS250kSimulation.R 1 300 3 6 1 0.3 0.4 31&
Rscript neighborGWAS250kSimulation.R 1 300 3 6 1 0.1 0.4 32&

Rscript neighborGWAS250kSimulation.R 0.25 300 3 6 1 0.6 0.8 33&
Rscript neighborGWAS250kSimulation.R 0.25 300 3 6 1 0.3 0.8 34&
Rscript neighborGWAS250kSimulation.R 0.25 300 3 6 1 0.3 0.4 35&
Rscript neighborGWAS250kSimulation.R 0.25 300 3 6 1 0.1 0.4 36&

#run at 50 SNPs
Rscript neighborGWAS250kSimulation.R 3 50 6 3 1 0.6 0.8 37&
Rscript neighborGWAS250kSimulation.R 3 50 6 3 1 0.3 0.8 38&
Rscript neighborGWAS250kSimulation.R 3 50 6 3 1 0.3 0.4 39&
Rscript neighborGWAS250kSimulation.R 3 50 6 3 1 0.1 0.4 40&

Rscript neighborGWAS250kSimulation.R 1 50 6 3 1 0.6 0.8 41&
Rscript neighborGWAS250kSimulation.R 1 50 6 3 1 0.3 0.8 42&
Rscript neighborGWAS250kSimulation.R 1 50 6 3 1 0.3 0.4 43&
Rscript neighborGWAS250kSimulation.R 1 50 6 3 1 0.1 0.4 44&

Rscript neighborGWAS250kSimulation.R 0.25 50 6 3 1 0.6 0.8 45&
Rscript neighborGWAS250kSimulation.R 0.25 50 6 3 1 0.3 0.8 46&
Rscript neighborGWAS250kSimulation.R 0.25 50 6 3 1 0.3 0.4 47&
Rscript neighborGWAS250kSimulation.R 0.25 50 6 3 1 0.1 0.4 48&

Rscript neighborGWAS250kSimulation.R 3 50 4 4 1 0.6 0.8 49&
Rscript neighborGWAS250kSimulation.R 3 50 4 4 1 0.3 0.8 50&
Rscript neighborGWAS250kSimulation.R 3 50 4 4 1 0.3 0.4 51&
Rscript neighborGWAS250kSimulation.R 3 50 4 4 1 0.1 0.4 52&

Rscript neighborGWAS250kSimulation.R 1 50 4 4 1 0.6 0.8 53&
Rscript neighborGWAS250kSimulation.R 1 50 4 4 1 0.3 0.8 54&
wait
Rscript neighborGWAS250kSimulation.R 1 50 4 4 1 0.3 0.4 55&
Rscript neighborGWAS250kSimulation.R 1 50 4 4 1 0.1 0.4 56&

Rscript neighborGWAS250kSimulation.R 0.25 50 4 4 1 0.6 0.8 57&
Rscript neighborGWAS250kSimulation.R 0.25 50 4 4 1 0.3 0.8 58&
Rscript neighborGWAS250kSimulation.R 0.25 50 4 4 1 0.3 0.4 59&
Rscript neighborGWAS250kSimulation.R 0.25 50 4 4 1 0.1 0.4 60&

Rscript neighborGWAS250kSimulation.R 3 50 3 6 1 0.6 0.8 61&
Rscript neighborGWAS250kSimulation.R 3 50 3 6 1 0.3 0.8 62&
Rscript neighborGWAS250kSimulation.R 3 50 3 6 1 0.3 0.4 63&
Rscript neighborGWAS250kSimulation.R 3 50 3 6 1 0.1 0.4 64&

Rscript neighborGWAS250kSimulation.R 1 50 3 6 1 0.6 0.8 65&
Rscript neighborGWAS250kSimulation.R 1 50 3 6 1 0.3 0.8 66&
Rscript neighborGWAS250kSimulation.R 1 50 3 6 1 0.3 0.4 67&
Rscript neighborGWAS250kSimulation.R 1 50 3 6 1 0.1 0.4 68&

Rscript neighborGWAS250kSimulation.R 0.25 50 3 6 1 0.6 0.8 69&
Rscript neighborGWAS250kSimulation.R 0.25 50 3 6 1 0.3 0.8 70&
Rscript neighborGWAS250kSimulation.R 0.25 50 3 6 1 0.3 0.4 71&
Rscript neighborGWAS250kSimulation.R 0.25 50 3 6 1 0.1 0.4 72&

#run at 10 SNPs
Rscript neighborGWAS250kSimulation.R 3 10 6 3 1 0.6 0.8 73&
Rscript neighborGWAS250kSimulation.R 3 10 6 3 1 0.3 0.8 74&
Rscript neighborGWAS250kSimulation.R 3 10 6 3 1 0.3 0.4 75&
Rscript neighborGWAS250kSimulation.R 3 10 6 3 1 0.1 0.4 76&

Rscript neighborGWAS250kSimulation.R 1 10 6 3 1 0.6 0.8 77&
Rscript neighborGWAS250kSimulation.R 1 10 6 3 1 0.3 0.8 78&
Rscript neighborGWAS250kSimulation.R 1 10 6 3 1 0.3 0.4 79&
Rscript neighborGWAS250kSimulation.R 1 10 6 3 1 0.1 0.4 80&

Rscript neighborGWAS250kSimulation.R 0.25 10 6 3 1 0.6 0.8 81&
wait
Rscript neighborGWAS250kSimulation.R 0.25 10 6 3 1 0.3 0.8 82&
Rscript neighborGWAS250kSimulation.R 0.25 10 6 3 1 0.3 0.4 83&
Rscript neighborGWAS250kSimulation.R 0.25 10 6 3 1 0.1 0.4 84&

Rscript neighborGWAS250kSimulation.R 3 10 4 4 1 0.6 0.8 85&
Rscript neighborGWAS250kSimulation.R 3 10 4 4 1 0.3 0.8 86&
Rscript neighborGWAS250kSimulation.R 3 10 4 4 1 0.3 0.4 87&
Rscript neighborGWAS250kSimulation.R 3 10 4 4 1 0.1 0.4 88&

Rscript neighborGWAS250kSimulation.R 1 10 4 4 1 0.6 0.8 89&
Rscript neighborGWAS250kSimulation.R 1 10 4 4 1 0.3 0.8 90&
Rscript neighborGWAS250kSimulation.R 1 10 4 4 1 0.3 0.4 91&
Rscript neighborGWAS250kSimulation.R 1 10 4 4 1 0.1 0.4 92&

Rscript neighborGWAS250kSimulation.R 0.25 10 4 4 1 0.6 0.8 93&
Rscript neighborGWAS250kSimulation.R 0.25 10 4 4 1 0.3 0.8 94&
Rscript neighborGWAS250kSimulation.R 0.25 10 4 4 1 0.3 0.4 95&
Rscript neighborGWAS250kSimulation.R 0.25 10 4 4 1 0.1 0.4 96&

Rscript neighborGWAS250kSimulation.R 3 10 3 6 1 0.6 0.8 97&
Rscript neighborGWAS250kSimulation.R 3 10 3 6 1 0.3 0.8 98&
Rscript neighborGWAS250kSimulation.R 3 10 3 6 1 0.3 0.4 99&
Rscript neighborGWAS250kSimulation.R 3 10 3 6 1 0.1 0.4 100&

Rscript neighborGWAS250kSimulation.R 1 10 3 6 1 0.6 0.8 101&
Rscript neighborGWAS250kSimulation.R 1 10 3 6 1 0.3 0.8 102&
Rscript neighborGWAS250kSimulation.R 1 10 3 6 1 0.3 0.4 103&
Rscript neighborGWAS250kSimulation.R 1 10 3 6 1 0.1 0.4 104&

Rscript neighborGWAS250kSimulation.R 0.25 10 3 6 1 0.6 0.8 105&
Rscript neighborGWAS250kSimulation.R 0.25 10 3 6 1 0.3 0.8 106&
Rscript neighborGWAS250kSimulation.R 0.25 10 3 6 1 0.3 0.4 107&
Rscript neighborGWAS250kSimulation.R 0.25 10 3 6 1 0.1 0.4 108&
wait

