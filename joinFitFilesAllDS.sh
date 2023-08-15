
#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -N snolab_monitoring
#$ -l h=compute-3-*.local
#$ -pe mpi 24 

source ~/env/damicm/bin/activate

# # run0d
# echo "Processing run0d"
# python joinFitsFiles.py -f /data/damic/snolab/upgrade/processed/science/run0d/avg/avg_skp_science_run0d_nrows210_ncols3100_nskips460_1_*_L.fits -o /data/damic/snolab/upgrade/processed/science/run0d/joined/1L
# python joinFitsFiles.py -f /data/damic/snolab/upgrade/processed/science/run0d/avg/avg_skp_science_run0d_nrows210_ncols3100_nskips460_1_*_U.fits -o /data/damic/snolab/upgrade/processed/science/run0d/joined/1U
# python joinFitsFiles.py -f /data/damic/snolab/upgrade/processed/science/run0d/avg/avg_skp_science_run0d_nrows210_ncols3100_nskips460_2_*_L.fits -o /data/damic/snolab/upgrade/processed/science/run0d/joined/2L
# python joinFitsFiles.py -f /data/damic/snolab/upgrade/processed/science/run0d/avg/avg_skp_science_run0d_nrows210_ncols3100_nskips460_2_*_U.fits -o /data/damic/snolab/upgrade/processed/science/run0d/joined/2U

# #run1a
# echo "Processing run1a"
# python joinFitsFiles.py -f /data/damic/snolab/upgrade/processed/science/run1a/avg/avg_skp_science_run1a_nrows210_ncols3100_nskips460_1_*_L.fits -o /data/damic/snolab/upgrade/processed/science/run1a/joined/1L
# python joinFitsFiles.py -f /data/damic/snolab/upgrade/processed/science/run1a/avg/avg_skp_science_run1a_nrows210_ncols3100_nskips460_1_*_U.fits -o /data/damic/snolab/upgrade/processed/science/run1a/joined/1U
# python joinFitsFiles.py -f /data/damic/snolab/upgrade/processed/science/run1a/avg/avg_skp_science_run1a_nrows210_ncols3100_nskips460_2_*_L.fits -o /data/damic/snolab/upgrade/processed/science/run1a/joined/2L
# python joinFitsFiles.py -f /data/damic/snolab/upgrade/processed/science/run1a/avg/avg_skp_science_run1a_nrows210_ncols3100_nskips460_2_*_U.fits -o /data/damic/snolab/upgrade/processed/science/run1a/joined/2U

# # run2b
echo "Processing run2b"
python joinFitsFiles.py -f /data/damic/snolab/upgrade/processed/science/run2b/avg/avg_skp_science_run2a_nrows210_ncols3300_nskips460_1_*_L.fits -o /data/damic/snolab/upgrade/processed/science/run2b/joined/1L
python joinFitsFiles.py -f /data/damic/snolab/upgrade/processed/science/run2b/avg/avg_skp_science_run2a_nrows210_ncols3300_nskips460_1_*_U.fits -o /data/damic/snolab/upgrade/processed/science/run2b/joined/1U
python joinFitsFiles.py -f /data/damic/snolab/upgrade/processed/science/run2b/avg/avg_skp_science_run2a_nrows210_ncols3300_nskips460_2_*_L.fits -o /data/damic/snolab/upgrade/processed/science/run2b/joined/2L
python joinFitsFiles.py -f /data/damic/snolab/upgrade/processed/science/run2b/avg/avg_skp_science_run2a_nrows210_ncols3300_nskips460_2_*_U.fits -o /data/damic/snolab/upgrade/processed/science/run2b/joined/2U


# # run3b
# echo "Processing run3b"
# python joinFitsFiles.py -f /data/damic/snolab/upgrade/processed/science/run3b/avg/avg_skp_science_run3b_nrows210_ncols3300_nskips460_1_*_L.fits -o /data/damic/snolab/upgrade/processed/science/run3b/joined/1L
# python joinFitsFiles.py -f /data/damic/snolab/upgrade/processed/science/run3b/avg/avg_skp_science_run3b_nrows210_ncols3300_nskips460_1_*_U.fits -o /data/damic/snolab/upgrade/processed/science/run3b/joined/1U
# python joinFitsFiles.py -f /data/damic/snolab/upgrade/processed/science/run3b/avg/avg_skp_science_run3b_nrows210_ncols3300_nskips460_2_*_L.fits -o /data/damic/snolab/upgrade/processed/science/run3b/joined/2L
# python joinFitsFiles.py -f /data/damic/snolab/upgrade/processed/science/run3b/avg/avg_skp_science_run3b_nrows210_ncols3300_nskips460_2_*_U.fits -o /data/damic/snolab/upgrade/processed/science/run3b/joined/2U


# # run4a
# echo "Processing run4a"
# python joinFitsFiles.py -f /data/damic/snolab/upgrade/processed/science/run4a/avg/avg_skp_science_run4a_nrows210_ncols3300_nskips460_1_*_L.fits -o /data/damic/snolab/upgrade/processed/science/run4a/joined/1L
# python joinFitsFiles.py -f /data/damic/snolab/upgrade/processed/science/run4a/avg/avg_skp_science_run4a_nrows210_ncols3300_nskips460_1_*_U.fits -o /data/damic/snolab/upgrade/processed/science/run4a/joined/1U
# python joinFitsFiles.py -f /data/damic/snolab/upgrade/processed/science/run4a/avg/avg_skp_science_run4a_nrows210_ncols3300_nskips460_2_*_L.fits -o /data/damic/snolab/upgrade/processed/science/run4a/joined/2L
# python joinFitsFiles.py -f /data/damic/snolab/upgrade/processed/science/run4a/avg/avg_skp_science_run4a_nrows210_ncols3300_nskips460_2_*_U.fits -o /data/damic/snolab/upgrade/processed/science/run4a/joined/2U

# echo "Processing run5a"
# python joinFitsFiles.py -f /data/damic/snolab/upgrade/processed/science/run5a/avg/avg_skp_science_run5a_nrows210_ncols3300_nskips460_1_*_L.fits -o /data/damic/snolab/upgrade/processed/science/run5a/joined/1L
# python joinFitsFiles.py -f /data/damic/snolab/upgrade/processed/science/run5a/avg/avg_skp_science_run5a_nrows210_ncols3300_nskips460_1_*_U.fits -o /data/damic/snolab/upgrade/processed/science/run5a/joined/1U
# python joinFitsFiles.py -f /data/damic/snolab/upgrade/processed/science/run5a/avg/avg_skp_science_run5a_nrows210_ncols3300_nskips460_2_*_L.fits -o /data/damic/snolab/upgrade/processed/science/run5a/joined/2L
# python joinFitsFiles.py -f /data/damic/snolab/upgrade/processed/science/run5a/avg/avg_skp_science_run5a_nrows210_ncols3300_nskips460_2_*_U.fits -o /data/damic/snolab/upgrade/processed/science/run5a/joined/2U

# run6a
#echo "Processing run5a"
#python joinFitsFiles.py -f /data/damic/snolab/upgrade/processed/science/run6a/avg/avg_skp_science_run6a_nrows210_ncols3300_nskips460_1_*_L.fits -o /data/damic/snolab/upgrade/processed/science/run6a/joined/1L
#python joinFitsFiles.py -f /data/damic/snolab/upgrade/processed/science/run6a/avg/avg_skp_science_run6a_nrows210_ncols3300_nskips460_1_*_U.fits -o /data/damic/snolab/upgrade/processed/science/run6a/joined/1U
#python joinFitsFiles.py -f /data/damic/snolab/upgrade/processed/science/run6a/avg/avg_skp_science_run6a_nrows210_ncols3300_nskips460_2_*_L.fits -o /data/damic/snolab/upgrade/processed/science/run6a/joined/2L
#python joinFitsFiles.py -f /data/damic/snolab/upgrade/processed/science/run6a/avg/avg_skp_science_run6a_nrows210_ncols3300_nskips460_2_*_U.fits -o /data/damic/snolab/upgrade/processed/science/run6a/joined/2U




