#################################################################################################################################
# nnSVG performances on simulations using increasing number of genes and cells. See Figure 2C,2D, and Supplementary Figure 3
#################################################################################################################################

# Genes
echo '2000*2000'
Rscript nnSVG.r /home/wangjue/Time_Mem_Cost_Test_v2/usage/Gene2000_Cell2000.csv > grid_nnSVG.txt
echo '2000*4000'
Rscript nnSVG.r /home/wangjue/Time_Mem_Cost_Test_v2/usage/Gene4000_Cell2000.csv >> grid_nnSVG.txt
echo '2000*6000'
Rscript nnSVG.r /home/wangjue/Time_Mem_Cost_Test_v2/usage/Gene6000_Cell2000.csv >> grid_nnSVG.txt
echo '2000*8000'
Rscript nnSVG.r /home/wangjue/Time_Mem_Cost_Test_v2/usage/Gene8000_Cell2000.csv >> grid_nnSVG.txt
echo '2000*10000'
Rscript nnSVG.r /home/wangjue/Time_Mem_Cost_Test_v2/usage/Gene10000_Cell2000.csv >> grid_nnSVG.txt
echo '2000*12000'
Rscript nnSVG.r /home/wangjue/Time_Mem_Cost_Test_v2/usage/Gene12000_Cell2000.csv >> grid_nnSVG.txt
echo '2000*14000'
Rscript nnSVG.r /home/wangjue/Time_Mem_Cost_Test_v2/usage/Gene14000_Cell2000.csv >> grid_nnSVG.txt
echo '2000*16000'
Rscript nnSVG.r /home/wangjue/Time_Mem_Cost_Test_v2/usage/Gene16000_Cell2000.csv >> grid_nnSVG.txt
echo '2000*18000'
Rscript nnSVG.r /home/wangjue/Time_Mem_Cost_Test_v2/usage/Gene18000_Cell2000.csv >> grid_nnSVG.txt
echo '2000*20000'
Rscript nnSVG.r /home/wangjue/Time_Mem_Cost_Test_v2/usage/Gene20000_Cell2000.csv >> grid_nnSVG.txt

#Cells
echo '1000*10000'
Rscript nnSVG.r /home/wangjue/Time_Mem_Cost_Test_v2/usage/Gene10000_Cell1000.csv >> grid_nnSVG.txt
echo '2000*10000'
Rscript nnSVG.r /home/wangjue/Time_Mem_Cost_Test_v2/usage/Gene10000_Cell2000.csv >> grid_nnSVG.txt
echo '3000*10000'
Rscript nnSVG.r /home/wangjue/Time_Mem_Cost_Test_v2/usage/Gene10000_Cell3000.csv >> grid_nnSVG.txt
echo '4000*10000'
Rscript nnSVG.r /home/wangjue/Time_Mem_Cost_Test_v2/usage/Gene10000_Cell4000.csv >> grid_nnSVG.txt
echo '5000*10000'
Rscript nnSVG.r /home/wangjue/Time_Mem_Cost_Test_v2/usage/Gene10000_Cell5000.csv >> grid_nnSVG.txt
echo '6000*10000'
Rscript nnSVG.r /home/wangjue/Time_Mem_Cost_Test_v2/usage/Gene10000_Cell6000.csv >> grid_nnSVG.txt
echo '7000*10000'
Rscript nnSVG.r /home/wangjue/Time_Mem_Cost_Test_v2/usage/Gene10000_Cell7000.csv >> grid_nnSVG.txt
echo '8000*10000'
Rscript nnSVG.r /home/wangjue/Time_Mem_Cost_Test_v2/usage/Gene10000_Cell8000.csv >> grid_nnSVG.txt
echo '9000*10000'
Rscript nnSVG.r /home/wangjue/Time_Mem_Cost_Test_v2/usage/Gene10000_Cell9000.csv >> grid_nnSVG.txt
echo '10000*10000'
Rscript nnSVG.r /home/wangjue/Time_Mem_Cost_Test_v2/usage/Gene10000_Cell10000.csv >> grid_nnSVG.txt




