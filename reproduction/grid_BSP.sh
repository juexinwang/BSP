#################################################################################################################################
# BSP performances on simulations using increasing number of genes and cells. See Figure 2C,2D, and Supplementary Figure 3
#################################################################################################################################

# Genes
echo '2000*2000'
python BSP.py --debugperform --oneFileTag --datasetName By_Genes --inputFile Gene2000_Cell2000.csv --inputDir /home/wangjue/Time_Mem_Cost_Test/ --outputDir /home/wangjue/svgGranularity/tmp/ > grid.txt
echo '2000*4000'
python BSP.py --debugperform --oneFileTag --datasetName By_Genes --inputFile Gene4000_Cell2000.csv --inputDir /home/wangjue/Time_Mem_Cost_Test/ --outputDir /home/wangjue/svgGranularity/tmp/ >> grid.txt
echo '2000*6000'
python BSP.py --debugperform --oneFileTag --datasetName By_Genes --inputFile Gene6000_Cell2000.csv --inputDir /home/wangjue/Time_Mem_Cost_Test/ --outputDir /home/wangjue/svgGranularity/tmp/ >> grid.txt
echo '2000*8000'
python BSP.py --debugperform --oneFileTag --datasetName By_Genes --inputFile Gene8000_Cell2000.csv --inputDir /home/wangjue/Time_Mem_Cost_Test/ --outputDir /home/wangjue/svgGranularity/tmp/ >> grid.txt
echo '2000*10000'
python BSP.py --debugperform --oneFileTag --datasetName By_Genes --inputFile Gene10000_Cell2000.csv --inputDir /home/wangjue/Time_Mem_Cost_Test/ --outputDir /home/wangjue/svgGranularity/tmp/ >> grid.txt
echo '2000*12000'
python BSP.py --debugperform --oneFileTag --datasetName By_Genes --inputFile Gene12000_Cell2000.csv --inputDir /home/wangjue/Time_Mem_Cost_Test/ --outputDir /home/wangjue/svgGranularity/tmp/ >> grid.txt
echo '2000*14000'
python BSP.py --debugperform --oneFileTag --datasetName By_Genes --inputFile Gene14000_Cell2000.csv --inputDir /home/wangjue/Time_Mem_Cost_Test/ --outputDir /home/wangjue/svgGranularity/tmp/ >> grid.txt
echo '2000*16000'
python BSP.py --debugperform --oneFileTag --datasetName By_Genes --inputFile Gene16000_Cell2000.csv --inputDir /home/wangjue/Time_Mem_Cost_Test/ --outputDir /home/wangjue/svgGranularity/tmp/ >> grid.txt
echo '2000*18000'
python BSP.py --debugperform --oneFileTag --datasetName By_Genes --inputFile Gene18000_Cell2000.csv --inputDir /home/wangjue/Time_Mem_Cost_Test/ --outputDir /home/wangjue/svgGranularity/tmp/ >> grid.txt
echo '2000*20000'
python BSP.py --debugperform --oneFileTag --datasetName By_Genes --inputFile Gene20000_Cell2000.csv --inputDir /home/wangjue/Time_Mem_Cost_Test/ --outputDir /home/wangjue/svgGranularity/tmp/ >> grid.txt

#Cells
echo '1000*10000'
python BSP.py --debugperform --oneFileTag --datasetName By_Cells --inputFile Gene10000_Cell1000.csv --inputDir /home/wangjue/Time_Mem_Cost_Test/ --outputDir /home/wangjue/svgGranularity/tmp/ >> grid.txt
echo '2000*10000'
python BSP.py --debugperform --oneFileTag --datasetName By_Cells --inputFile Gene10000_Cell2000.csv --inputDir /home/wangjue/Time_Mem_Cost_Test/ --outputDir /home/wangjue/svgGranularity/tmp/ >> grid.txt
echo '3000*10000'
python BSP.py --debugperform --oneFileTag --datasetName By_Cells --inputFile Gene10000_Cell3000.csv --inputDir /home/wangjue/Time_Mem_Cost_Test/ --outputDir /home/wangjue/svgGranularity/tmp/ >> grid.txt
echo '4000*10000'
python BSP.py --debugperform --oneFileTag --datasetName By_Cells --inputFile Gene10000_Cell4000.csv --inputDir /home/wangjue/Time_Mem_Cost_Test/ --outputDir /home/wangjue/svgGranularity/tmp/ >> grid.txt
echo '5000*10000'
python BSP.py --debugperform --oneFileTag --datasetName By_Cells --inputFile Gene10000_Cell5000.csv --inputDir /home/wangjue/Time_Mem_Cost_Test/ --outputDir /home/wangjue/svgGranularity/tmp/ >> grid.txt
echo '6000*10000'
python BSP.py --debugperform --oneFileTag --datasetName By_Cells --inputFile Gene10000_Cell6000.csv --inputDir /home/wangjue/Time_Mem_Cost_Test/ --outputDir /home/wangjue/svgGranularity/tmp/ >> grid.txt
echo '7000*10000'
python BSP.py --debugperform --oneFileTag --datasetName By_Cells --inputFile Gene10000_Cell7000.csv --inputDir /home/wangjue/Time_Mem_Cost_Test/ --outputDir /home/wangjue/svgGranularity/tmp/ >> grid.txt
echo '8000*10000'
python BSP.py --debugperform --oneFileTag --datasetName By_Cells --inputFile Gene10000_Cell8000.csv --inputDir /home/wangjue/Time_Mem_Cost_Test/ --outputDir /home/wangjue/svgGranularity/tmp/ >> grid.txt
echo '9000*10000'
python BSP.py --debugperform --oneFileTag --datasetName By_Cells --inputFile Gene10000_Cell9000.csv --inputDir /home/wangjue/Time_Mem_Cost_Test/ --outputDir /home/wangjue/svgGranularity/tmp/ >> grid.txt
echo '10000*10000'
python BSP.py --debugperform --oneFileTag --datasetName By_Cells --inputFile Gene10000_Cell10000.csv --inputDir /home/wangjue/Time_Mem_Cost_Test/ --outputDir /home/wangjue/svgGranularity/tmp/ >> grid.txt




