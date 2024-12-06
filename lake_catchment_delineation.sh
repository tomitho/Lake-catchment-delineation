#!/bin/bash

#SBATCH --job-name=lakejob
#SBATCH --output /gpfs/gibbs/project/sbsc/tt572/lakes/monitor/slurm-%A_%a-%N.out
#SBATCH --partition=scavenge
#SBATCH --time 00:30:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu 15G
#SBATCH --requeue
#SBATCH --array=1-5754

module purge
module load PKTOOLS/2.6.7.6-foss-2022b
module load GRASS/8.2.0-foss-2022b
module load GDAL/3.6.2-foss-2022b
module load GnuTLS/3.7.8-GCCcore-12.2.0
module load PostgreSQL/15.2-GCCcore-12.2.0


# Specify the path to the config file, contains all lake ids being processed in the array (max. 10000 per array)
config=/gpfs/gibbs/project/sbsc/tt572/lakes/config_files/config.txt

# Extract the lake ID for the current $SLURM_ARRAY_TASK_ID
LK=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)

export LK=$LK

# global file of basins
export GLBASINS=/gpfs/gibbs/pi/hydro/hydro/dataproces/MERIT_HYDRO/lbasin_tiles_final20d_ovr/all_lbasin.tif
# global basin of computational units
export GLCOMPUNITS=/gpfs/gibbs/pi/hydro/hydro/dataproces/MERIT_HYDRO/lbasin_compUnit_overview/lbasin_compUnit.tif
# path  - basins in each computational unit
export CUBASINS=/gpfs/gibbs/pi/hydro/hydro/dataproces/MERIT_HYDRO/CompUnit_lbasin
# global vrt file of streams
export STREAM=/home/tt572/vrts/all_streams_uniqID.vrt
# global vrt file of flow accumulation
export FLOW=/home/tt572/vrts/flow_global.vrt
# global vrt file of direction
export DIRECTION=/gpfs/gibbs/pi/hydro/hydro/dataproces/MERIT_HYDRO/dir_tiles_final20d_1p/all_dir_dis.vrt
# path where output data should be directed 
export DIR=/gpfs/gibbs/project/sbsc/tt572/lakes
# create for each computing node a tmp folder to not clock the tmp directory when large amounts of arraytask accessing tmp/
mkdir /gpfs/gibbs/project/sbsc/tt572/lakes/tmp/$(hostname)/
mkdir /gpfs/gibbs/project/sbsc/tt572/lakes/tmp/$(hostname)/$SLURM_JOB_ID
# path to temporal /compilation
export tmp=/gpfs/gibbs/project/sbsc/tt572/lakes/tmp/$(hostname)/$SLURM_JOB_ID
# path to regional unit folder
export RE=/gpfs/gibbs/project/sbsc/tt572/lakes/out
# path to regional unit look up
export LO=/gpfs/gibbs/project/sbsc/tt572/lakes/lake_area.txt
# software path for MSPA analysis 
export GWB=/gpfs/gibbs/project/sbsc/tt572/software/GWB1.9.4/GWB

###############################################################################
###############################################################################
###############################################################################

echo "START THE SCRIPT"

LAK=$(awk -v lk="$LK" '$1 == lk {print $2}' $LO)

### print regional unit ###

echo --------------------------------------------
echo "Lake ${LK} falls into Regional Unit ${LAK}"
echo --------------------------------------------


# select lake shapefile
    ogr2ogr $tmp/lake_${LK}.shp \
	     $DIR/lakes.gpkg \
       -sql "SELECT * FROM HydroLAKES_polys_v10 WHERE Hylak_id = ${LK}"

        export LAKE=$tmp/lake_${LK}.shp

    # check basename of polygon file
    export pn=$(basename $LAKE .shp)


echo
echo ------------------------
echo "Run Script for Lake ${LK}"
echo ------------------------
echo

bd=$(awk -v lk="$LK" '$1 == lk {print $4}' $LO)

    ## calculate buffer
    ogr2ogr -dialect sqlite -sql "SELECT ST_Buffer(Geometry, ${bd}) FROM ${pn}" \
        $tmp/buffer_${LK}.shp \
        $tmp/lake_${LK}.shp

##  lake file (e.g. shp)
export LAKE=$tmp/buffer_${LK}.shp

# create folders for this lake based on id
if ! [[ -d $DIR/out/lake_${LK}/supplementary ]]
then
mkdir -p $DIR/out/lake_${LK}
fi
export OUTDIR=$DIR/out/lake_${LK}

# supplementary data
if ! [[ -d $DIR/out/lake_${LK}/supplementary ]]
then
mkdir $DIR/out/lake_${LK}/supplementary
fi
export SUPPL=$DIR/out/lake_${LK}/supplementary

# check basename of polygon file
export pn=$(basename $LAKE .shp)

if ! [[ -f $DIR/out/lake_${Lk}/supplementary/lake_${LK}.tif ]]
then
# check bbox and extent to integer number
# to secure spatial match with Hydrography90m dataset
EXTENSION=($( ogrinfo  $LAKE -so -al | grep Extent \
    | grep -Eo '[+-]?[0-9]+([.][0-9]+)?' ))

temp=${EXTENSION[0]}
if (($(bc <<< "$temp < 0")))
then
    xmin=$(echo $temp | awk '{print int($1)-1}')
else
    xmin=$(echo $temp | awk '{print int($1)}')
fi

temp=${EXTENSION[2]}
if (($(bc <<< "$temp < 0")))
then
    xmax=$(echo $temp | awk '{print int($1)}')
else
    xmax=$(echo $temp | awk '{print int($1)+1}')
fi

temp=${EXTENSION[1]}
if (($(bc <<< "$temp < 0")))
then
    ymin=$(echo $temp | awk '{print int($1)-1}')
else
    ymin=$(echo $temp | awk '{print int($1)}')
fi

temp=${EXTENSION[3]}
if (($(bc <<< "$temp < 0")))
then
    ymax=$(echo $temp | awk '{print int($1)}')
else
    ymax=$(echo $temp | awk '{print int($1)+1}')
fi


## Rasterize lake to use later as the mask layer

# add column as reference for raster creation
ogrinfo $LAKE -sql "ALTER TABLE $pn  ADD COLUMN diss INTEGER"
ogrinfo $LAKE -dialect SQLite -sql "UPDATE $pn SET diss = 1"

gdal_rasterize -a_srs EPSG:4326  -at -a diss -l $pn \
    -tr 0.000833333333333 -0.000833333333333 \
    -te $xmin $ymin $xmax $ymax -a_nodata 0  \
    -co COMPRESS=DEFLATE -co ZLEVEL=9 -ot Byte \
    $LAKE $tmp/lake_${LK}cp.tif

echo
echo ------------------------------------------------------
echo "EXTRACTING THE COMPUNIT AND BASINS IDS FOR LAKE ${LK}"
echo ------------------------------------------------------
echo

# crop by removing NA data
grass -f --text --tmp-location $tmp/lake_${LK}cp.tif <<'EOF'
    r.external -o input=$tmp/lake_${LK}cp.tif  output=out
    g.region zoom=out
    r.out.gdal -cm input=out out=$tmp/lake_${LK}rm.tif \
    format=GTiff type=Byte createopt="COMPRESS=DEFLATE,BIGTIFF=YES" nodata=0 --overwrite
EOF

# add again a couple of pixels all around
# to avoid border breaks when running MSPA
xmin=$(pkinfo -i $tmp/lake_${LK}rm.tif -te | awk '{print $2-0.001666667}')
ymin=$(pkinfo -i $tmp/lake_${LK}rm.tif -te | awk '{print $3-0.001666667}')
xmax=$(pkinfo -i $tmp/lake_${LK}rm.tif -te | awk '{print $4+0.001666667}')
ymax=$(pkinfo -i $tmp/lake_${LK}rm.tif -te | awk '{print $5+0.001666667}')

gdalwarp -te $xmin $ymin $xmax $ymax -tr 0.000833333333333 -0.000833333333333 \
    -co COMPRESS=DEFLATE -co ZLEVEL=9 -ot Byte \
    $tmp/lake_${LK}rm.tif $SUPPL/lake_${LK}.tif


# crop basin raster file to the lake extention

gdal_translate \
    -projwin $( pkinfo -i $SUPPL/lake_${LK}.tif -bb | grep -Eo '[+-]?[0-9]+([.][0-9]+)?' ) \
    -co COMPRESS=LZW -co ZLEVEL=9 $GLBASINS $tmp/extentLB_${LK}.tif

# crop compunit raster file to the lake extention
gdal_translate \
    -projwin $( pkinfo -i $SUPPL/lake_${LK}.tif -bb | grep -Eo '[+-]?[0-9]+([.][0-9]+)?' ) \
    -co COMPRESS=LZW -co ZLEVEL=9 $GLCOMPUNITS $tmp/extentCU_${LK}.tif

# mask the basins raster to the lake boundaries to identify the basin IDs
# that fully overlap with lake
pksetmask -i $tmp/extentLB_${LK}.tif \
    -m $SUPPL/lake_${LK}.tif -msknodata 0 -nodata 0 \
    -o $tmp/maskLB_${LK}.tif -co COMPRESS=LZW -co ZLEVEL=9

# identify basin IDs
echo "$LK $(pkstat -i $tmp/maskLB_${LK}.tif -hist \
    | awk '$2 > 0 && $1 > 0 {print $1}' | sed -z 's/\n/ /g')" \
    > $tmp/ref_lbasinIDs_${LK}.txt

# mask the compunit raster to the lake boundaries to identify the compUnit IDs
# that fully overlap with lake
pksetmask -i $tmp/extentCU_${LK}.tif \
    -m $SUPPL/lake_${LK}.tif -msknodata 0 -nodata 0 \
    -o $tmp/maskCU_${LK}.tif -co COMPRESS=LZW -co ZLEVEL=9

# identify CompUnit IDs
echo "$LK $(pkstat -i $tmp/maskCU_${LK}.tif -hist \
    | awk '$2 > 0 && $1 > 0 {print $1}' | sed -z 's/\n/ /g')" \
    > $tmp/ref_CompUnitIDs_${LK}.txt

# remove tmp files
rm $tmp/maskCU_${LK}.tif $tmp/maskLB_${LK}.tif \
$tmp/extentCU_${LK}.tif $tmp/extentLB_${LK}.tif \
$tmp/lake_${LK}rm.tif $tmp/lake_${LK}cp.tif



###############################################################################
###############################################################################
###############################################################################

echo ------------------------------------------------------
echo "EXTRACTING EXTENT OF ALL BASINS AROUND THE LAKE ${LK}"
echo ------------------------------------------------------


# check if the lake overlaps with more than two CompUnits
if [[ $(awk '{print NF}' $tmp/ref_CompUnitIDs_${LK}.txt) -gt 2 ]];
then

echo ------------------------------------------------------
echo "Its a MULTI Lake ${LK}"
echo ------------------------------------------------------


    # get the IDs of the COmpUnits
    CUID=($(cut -d" " -f2- $tmp/ref_CompUnitIDs_${LK}.txt))

    # Join Computational Units as VRT
    # create file with list of files to join
    for i in ${CUID[@]}; do find $CUBASINS -name "lbasin_${i}_msk.tif"; done \
        > $tmp/my_list_${LK}.txt
    # create the vrt
    gdalbuildvrt $tmp/CompUnits_${LK}.vrt \
        -input_file_list $tmp/my_list_${LK}.txt

    rm $tmp/my_list_${LK}.txt

    # make list of lbasins IDs as an arrray
    export  MBID=( $(cut -d" " -f2- $tmp/ref_lbasinIDs_${LK}.txt) )

   # reclassify raster with lbasins of interest (1) else is (0)
    # identify max value (lbasin ID)
    MAX=$(gdalinfo -mm $tmp/CompUnits_${LK}.vrt | grep Computed \
        | awk -F, '{print $2}')
    # create classification table
    col1=($(seq 0 ${MAX%.*}))
    col2=($(printf '0%.0s\n' $(eval "echo {0.."$((${MAX%.*}))"}" )))
    for i in ${MBID[@]}; do col2[$i]=1; done
    paste -d " " \
        <(printf "%s\n" ${col1[@]:1:${#col1[@]}}) \
        <(printf "%s\n" ${col2[@]:1:${#col2[@]}}) \
        > $tmp/reclass_code_${LK}.txt

    # run reclassification
    pkreclass -i $tmp/CompUnits_${LK}.vrt \
        -o $tmp/lbasin_reclass_${LK}.vrt \
        --code $tmp/reclass_code_${LK}.txt

    rm $tmp/reclass_code_${LK}.txt


grass -f --text --tmp-location $tmp/lbasin_reclass_${LK}.vrt <<'EOF'
    r.external -o input=$tmp/lbasin_reclass_${LK}.vrt output=out
    g.region zoom=out
    g.region -p
    r.out.gdal input=out out=${SUPPL}/ALLbasins_${LK}.vrt \
    format=VRT type=Byte --overwrite
    region_info=$(g.region -w | awk -F '[=,]' '{print $2,$3,$4,$5}')
    echo "$LK $region_info" > $tmp/lakes_ref_extent_${LK}.txt
EOF

rm $SUPPL/ALLbasins_${LK}.tif
    
    REG=multi
    if ! [[ -d $OUTDIR/reg_unit_${REG} ]]
    then
    mkdir $RE/reg_unit_${REG}
    fi
    
    rm  $tmp/lbasin_reclass_${LK}.vrt  $tmp/CompUnits_${LK}.vrt

else

CUID=($(cut -d" " -f2- $tmp/ref_CompUnitIDs_${LK}.txt))


# identify the computational unit where the basin(s) are
compunit=$(awk '{print $2}' $tmp/ref_CompUnitIDs_${LK}.txt)
export compunit

mkdir $SUPPL/$SLURM_JOB_ID

    # Join Computational Units as VRT
cp $CUBASINS/lbasin_${compunit}_msk.tif $SUPPL/$SLURM_JOB_ID

    # make list of lbasins IDs as an arrray
export  MBID=( $(cut -d" " -f2- $tmp/ref_lbasinIDs_${LK}.txt) )

   # reclassify raster with lbasins of interest (1) else is (0)
    # identify max value (lbasin ID)
    MAX=$(gdalinfo -mm $SUPPL/$SLURM_JOB_ID/lbasin_${compunit}_msk.tif | grep Computed \
        | awk -F, '{print $2}')
    # create classification table
    col1=($(seq 0 ${MAX%.*}))
    col2=($(printf '0%.0s\n' $(eval "echo {0.."$((${MAX%.*}))"}" )))
    for i in ${MBID[@]}; do col2[$i]=1; done
    paste -d "=" \
        <(printf "%s\n" ${col1[@]:1:${#col1[@]}}) \
        <(printf "%s\n" ${col2[@]:1:${#col2[@]}}) \
        > $tmp/reclass_code_${LK}.txt

    

grass -f --text --tmp-location $SUPPL/$SLURM_JOB_ID/lbasin_${compunit}_msk.tif <<'EOF'
   r.in.gdal input=$SUPPL/$SLURM_JOB_ID/lbasin_${compunit}_msk.tif output=mask
   g.region rast=mask
   r.reclass input=mask output=out rules=$tmp/reclass_code_${LK}.txt --overwrite
   if g.list type=rast pattern=out | grep -q "^out$"; then
   r.out.gdal input=out output=$SUPPL/$SLURM_JOB_ID/lbasin_reclass_${LK}.tif \
   format=GTiff type=Int32 --overwrite -c
   else
        echo "ERROR: Reclassification did not produce an output raster."
    fi
EOF


echo $SUPPL/file.tif
echo $SUPPL/file.vrt


gdal_translate -of VRT $SUPPL/$SLURM_JOB_ID/lbasin_reclass_${LK}.tif $SUPPL/$SLURM_JOB_ID/lbasin_reclass_${LK}.vrt


rm $tmp/reclass_code_${LK}.txt
rm $SUPPL/$SLURM_JOB_ID/lbasin_${compunit}_msk.tif
rm $SUPPL/$SLURM_JOB_ID/lbasin_reclass_${LK}.tif


grass -f --text --tmp-location $SUPPL/$SLURM_JOB_ID/lbasin_reclass_${LK}.vrt <<'EOF'
   r.external -o input=$SUPPL/$SLURM_JOB_ID/lbasin_reclass_${LK}.vrt output=out
   g.region zoom=out
   g.region -p
   r.out.gdal input=out out=${SUPPL}/ALLbasins_${LK}.tif \
   format=GTiff type=Byte createopt="COMPRESS=DEFLATE,BIGTIFF=YES" --overwrite
   region_info=$(g.region -w | awk -F '[=,]' '{print $2,$3,$4,$5}')
   echo "$LK $region_info" > $tmp/lakes_ref_extent_${LK}.txt
EOF

rm $SUPPL/ALLbasins_${LK}.tif
    
    REG=${compunit}
    if ! [[ -d $OUTDIR/reg_unit_${REG} ]]
    then
    mkdir $RE/reg_unit_${REG}
    fi
    
    rm $SUPPL/$SLURM_JOB_ID/lbasin_reclass_${LK}.vrt

fi

# extention to consider for upstream creation
export EXT=$(awk '{print $2, $3, $4, $5}' $tmp/lakes_ref_extent_${LK}.txt)

else
echo "Nothing to see here"
fi
#rm $OUTDIR/*ref_*.txt

###############################################################################
###############################################################################
###############################################################################

echo ------------------------------------------------------
echo "Do the edge analysis ${LK}"
echo ------------------------------------------------------

    mv $tmp/lake_${LK}.* $DIR/out/lake_${LK}/supplementary
    mv $tmp/buffer_${LK}.* $DIR/out/lake_${LK}/supplementary

SUPPL=$DIR/out/lake_${LK}/supplementary
export SUPPL
mkdir $DIR/GWB/$(hostname)/
mkdir $DIR/GWB/$(hostname)/lake_${LK}
cp -r $GWB $DIR/GWB/$(hostname)/lake_${LK}

pkreclass -co COMPRESS=LZW -co ZLEVEL=9 -ot Byte -of GTiff \
        -i $SUPPL/lake_${LK}.tif \
        -o $DIR/GWB/$(hostname)/lake_${LK}/GWB/input/lake_${LK}.tif \
        -c 1 -r 2 -c 0 -r 1

pkreclass -co COMPRESS=LZW -co ZLEVEL=9 -ot Byte -of GTiff \
        -i $DIR/GWB/$(hostname)/lake_${LK}/GWB/input/lake_${LK}.tif \
        -o $DIR/GWB/$(hostname)/lake_${LK}/GWB/input/lake_${LK}.tif \
        -c 255 -r 1


cd $DIR/GWB/$(hostname)/lake_${LK}/GWB
./GWB_MSPA -i=$DIR/GWB/$(hostname)/lake_${LK}/GWB/input -o=$DIR/GWB/$(hostname)/lake_${LK}/GWB/output
cd $DIR

SUPPL=$DIR/out/lake_${LK}/supplementary
cp $DIR/GWB/$(hostname)/lake_${LK}/GWB/output/lake_${LK}_mspa/lake_${LK}_8_1_0_1.tif $SUPPL

rm -r $DIR/GWB/$(hostname)/lake_${LK}*

### reclassification of MSPA output to retain only the outer borders
pkreclass -co COMPRESS=LZW -co ZLEVEL=9 \
    -i $SUPPL/lake_${LK}_8_1_0_1.tif \
    -o $SUPPL/mspa_${LK}.tif --code $DIR/mspa_reclass_code.txt

# warp to avoid pixel missmatch
gdalwarp -tr 0.000833333333333 0.000833333333333 -tap $SUPPL/mspa_${LK}.tif $SUPPL/mspa_${LK}_l.tif -overwrite

# crop streams to the extension of interest
gdal_translate -co COMPRESS=LZW -co ZLEVEL=9 \
    -projwin $( pkinfo -i $SUPPL/mspa_${LK}_l.tif -bb | grep  -Eo '[+-]?[0-9]+([.][0-9]+)?' ) \
    -tr 0.000833333333333 0.000833333333333 $STREAM $SUPPL/stream_${LK}.tif

###  Multiply to identify intersections
gdal_calc.py -A $SUPPL/stream_${LK}.tif -B $SUPPL/mspa_${LK}_l.tif \
    --outfile=$tmp/overlapping_${LK}.tif --calc="A*B" \
    --NoDataValue=0 --overwrite

### crop flow to extension of interest
gdal_translate -co COMPRESS=LZW -co ZLEVEL=9 \
    -projwin $( pkinfo -i $SUPPL/mspa_${LK}_l.tif -bb | grep  -Eo '[+-]?[0-9]+([.][0-9]+)?' ) \
    $FLOW $SUPPL/flow_${LK}.tif

cp $SUPPL/mspa_${LK}_l.tif $SUPPL/mspa_${LK}.tif
rm $SUPPL/mspa_${LK}_l.tif

grass  -f --text --tmp-location $tmp/overlapping_${LK}.tif <<'EOF'
r.external -o input=$tmp/overlapping_${LK}.tif  output=out
r.external -o input=$SUPPL/flow_${LK}.tif output=flow
r.neighbors input=flow output=flow_max method=maximum
r.neighbors input=flow output=flow_ave method=average

r.mapcalc --o "flowInt = if(isnull(out), null() , flow)"
r.mapcalc --o "flowMax = if(isnull(out), null() , flow_max)"
r.mapcalc --o "flowAve = if(isnull(out), null() , flow_ave)"

r.out.xyz --overwrite input=out output=$tmp/coord_lake_${LK}.txt separator=space
r.out.xyz --overwrite input=flowInt output=$tmp/coord_flowInt_${LK}.txt separator=space
r.out.xyz --overwrite input=flowMax output=$tmp/coord_flowMax_${LK}.txt separator=space
r.out.xyz --overwrite input=flowAve output=$tmp/coord_flowAve_${LK}.txt separator=space
EOF

paste -d " " <(seq 1 $(wc -l < $tmp/coord_lake_${LK}.txt))  \
    $tmp/coord_lake_${LK}.txt \
    <(awk '{printf "%.3f\n", $3}' $tmp/coord_flowInt_${LK}.txt)   \
    <(awk '{printf "%.3f\n", $3}' $tmp/coord_flowMax_${LK}.txt)   \
    <(awk '{printf "%.3f\n", $3}' $tmp/coord_flowAve_${LK}.txt)   \
    > $OUTDIR/coord_lake_${LK}.txt

rm $tmp/coord_lake_${LK}.txt $tmp/coord_flowInt_${LK}.txt
rm $tmp/coord_flowMax_${LK}.txt $tmp/coord_flowAve_${LK}.txt
rm $tmp/overlapping_${LK}.tif


echo ---------------------
echo REMOVING DUPLICATES FOR LAKE ${LK}
echo ---------------------

# Identify streamIDs that have more than one duplicate
DUP=$(awk '{print $4}' $OUTDIR/coord_lake_${LK}.txt | sort | uniq -c | awk '$1 > 1 {print $2}')

# Loop to go through each of the stream IDs and identify the coordIDs of the duplicates with the lowest flow accumulation. The purpose is to remove those from the main table and leave only one streamID with the highest accumulation value (using as comparison the Maxumin flow $6)
###  Procedure to retain the intersection with the minimum accumulation value (for the maximum replace tail with head)
RM=$(for LINE in $DUP
do
awk -v line=$LINE '$4 == line {print $1, $6}' $OUTDIR/coord_lake_${LK}.txt \
    | LC_ALL=C sort -k2 -g \
    | head -n $(echo "$(awk -v line=$LINE '$4 == line' $OUTDIR/coord_lake_${LK}.txt | wc -l)" - 1 | bc ) \
    | awk '{print $1}'
done)

# Move through each coordID duplicate and remove it from the main table
for ROW in $RM
do
awk -i inplace -v row=$ROW '$1 == row {next} {print}' $OUTDIR/coord_lake_${LK}.txt
done

### sort coordinate table after flow accumulation value which was used to select outlet ###

sort -k7 -n -r $OUTDIR/coord_lake_${LK}.txt > $OUTDIR/coord_lake_${LK}_sort.txt
paste -d " " <(seq 1 $(wc -l < $OUTDIR/coord_lake_${LK}_sort.txt))  \
  $OUTDIR/coord_lake_${LK}_sort.txt \
  > $OUTDIR/coord_lake_${LK}_tmp.txt

awk -v var=$LK '$(NF+1) = var' $OUTDIR/coord_lake_${LK}_tmp.txt \
> $OUTDIR/coord_lake_${LK}_pro.txt

rm $OUTDIR/coord_lake_${LK}.txt

awk '{print $1,$9,$3,$4,$5,$6,$7,$8}' "$OUTDIR/coord_lake_${LK}_pro.txt" \
> $OUTDIR/coord_lake_${LK}.txt

rm $OUTDIR/coord_lake_${LK}_pro.txt
rm $OUTDIR/coord_lake_${LK}_tmp.txt

cp $OUTDIR/coord_lake_${LK}.txt $SUPPL

###############################################################################
###############################################################################
###############################################################################

echo --------------------------------
echo SELECT THE OUTLET FOR LAKE ${LK}
echo --------------------------------

export TOP=1

sort -su -t, -k7 $OUTDIR/coord_lake_${LK}.txt > $OUTDIR/coord_lake_${LK}_t.txt
mv  $OUTDIR/coord_lake_${LK}_t.txt $OUTDIR/coord_lake_${LK}_top${TOP}.txt

echo ---------------------------------
echo CREATING UPSTREAM BASINS IN INTERSECTIONS FOR LAKE ${LK}
echo ---------------------------------

# extention to consider for upstream creation
export EXT=$(awk '{print $2, $3, $4, $5}' $tmp/lakes_ref_extent_${LK}.txt)
export DATACOORD=$OUTDIR/coord_lake_${LK}_top${TOP}.txt

## create lake raster layer with whole extention (add NAs)
gdalwarp  -co COMPRESS=LZW -co ZLEVEL=9 -te $EXT \
    $SUPPL/lake_${LK}.tif \
    $tmp/lake_${LK}_NA.tif -overwrite

gdalwarp  -co COMPRESS=LZW -co ZLEVEL=9 -te $EXT \
    $DIRECTION $tmp/dir_${LK}_NA.tif -overwrite


grass  -f --text --tmp-location $tmp/lake_${LK}_NA.tif <<'EOF'

# read binary lake layer with same extent as direction
r.external --o input=$tmp/lake_${LK}_NA.tif  output=lake

#  read direction map
r.external --o -a input=$tmp/dir_${LK}_NA.tif output=dir

export COD=1

## calculate the sub-basin
r.water.outlet --overwrite input=dir output=bf_${COD} \
coordinates=$(cat $DATACOORD | awk -v coord=${COD} 'BEGIN{OFS=",";} $1==coord {print $3,$4}')

##  Export the basin as tif file
r.out.gdal --o -f -c -m  createopt="COMPRESS=DEFLATE,ZLEVEL=9" \
type=Byte  format=GTiff nodata=0 \
input=bf_${COD} output=$tmp/basin_lake_${LK}_coord_${COD}.tif
EOF

COD=$(cat $DATACOORD | awk '{print $1}' $DATACOORD)

# remove temporal files
rm $tmp/lake_${LK}_NA.tif $tmp/dir_${LK}_NA.tif

#-----------------------------

mkdir  $OUTDIR/upstreamBasins
mv $tmp/basin_lake_${LK}_coord_${COD}.tif $OUTDIR/upstreamBasins


### Convert table of intersections to csv
cat $OUTDIR/coord_lake_${LK}_top${TOP}.txt | awk 'NR == 1 { print "outlet_id lake_id XLon YLat streamID flowpix flowmax flowmean" }; 1' \
  > $OUTDIR/coord_lake_${LK}_top${TOP}2.txt
awk '{gsub(/ /,","); print}' $OUTDIR/coord_lake_${LK}_top${TOP}2.txt \
    > $OUTDIR/coord_lake_${LK}_top${TOP}.csv

# create gpkg of outlets
# create cvs table of intersections
ogr2ogr -f "GPKG" \
 -oo X_POSSIBLE_NAMES=XLon -oo Y_POSSIBLE_NAMES=YLat  -oo AUTODETECT_TYPE=YES \
 $SUPPL/outlets_${LK}.gpkg $OUTDIR/coord_lake_${LK}_top${TOP}.csv

rm $OUTDIR/coord_lake_${LK}_top${TOP}2.txt
rm $OUTDIR/coord_lake_${LK}_top${TOP}.csv


# clean

rm $OUTDIR/outlets_${LK}.csv
rm $OUTDIR/ids.txt
rm $tmp/Lake_${LK}_allupst.vrt
rm -r $tmp
rm $tmp/lakes_ref_extent_${LK}.txt
rm $tmp/ref_CompUnitIDs_${LK}.txt
rm $tmp/ref_lbasinIDs_${LK}.txt
rm $tmp/lbasin_${compunit}_msk.tif
rm $OUTDIR/coord_lake_${LK}_top1.txt
rm $OUTDIR/supplementary/ALLbasins_${LK}.tif
rm $OUTDIR/supplementary/flow_${LK}.tif
rm $OUTDIR/supplementary/lake_${LK}.dbf
rm $OUTDIR/supplementary/lake_${LK}.prj
rm $OUTDIR/supplementary/lake_${LK}.shp
rm $OUTDIR/supplementary/lake_${LK}.shx
rm $OUTDIR/supplementary/stream_${LK}.tif
rm $OUTDIR/supplementary/mspa_${LK}.tif
rm $OUTDIR/supplementary/lake_${LK}_8_1_0_1.tif
rm $OUTDIR/coord_lake_${LK}*
rm $OUTDIR/supplementary/buffer_${LK}*

# show where folder is saved
echo ---------------------------------------------------------------
echo LAKE IS ZIPPED AND SAVED UNDER "$RE/reg_unit_${REG}/lake_${LK}"
echo ---------------------------------------------------------------

zip -r -j $RE/reg_unit_${REG}/lake_${LK}.zip $OUTDIR/supplementary $OUTDIR/upstreamBasins
rm -r $OUTDIR

# done
echo "End da Script ${LK}"
echo

exit