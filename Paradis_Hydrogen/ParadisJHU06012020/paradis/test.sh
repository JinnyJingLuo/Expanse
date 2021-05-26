#! /bin/bash
dir_seq='FR-source2'
dir_vloop_arr=('100' '200' '300' '400')
# dir_shearloop_arr=('20')
for i in ${!dir_vloop_arr[@]}
do 
    if [ ! -d "Cases/${dir_vloop_arr[$i]}_FRsource" ]
    then
        mkdir -p "Cases/${dir_vloop_arr[$i]}_FRsource"
    #mkdir "$"Cell1$counter""
    fi
    cp frank_read_1.data "Cases/${dir_vloop_arr[$i]}FRsource/frank_read_${dir_vloop_arr[$i]}.data"
    # cp S_wall_data "Cases/${dir_vloop_arr[$i]}_${dir_shearloop_arr[$j]}FRsource"
    # cp S_interior_data "Cases/${dir_vloop_arr[$i]}_${dir_shearloop_arr[$j]}FRsource"
    cp data_FR_Jing2.py "Cases/${dir_vloop_arr[$i]}_FRsource"
    cp paradisJHU "Cases/${dir_vloop_arr[$i]}_FRsource"
    cp frank_read_1.ctrl "Cases/${dir_vloop_arr[$i]}_FRsource/frank_read_${dir_vloop_arr[$i]}.ctrl"
    cd "Cases/${dir_vloop_arr[$i]}_${dir_shearloop_arr[$j]}FRsource"
    sed -i "s/rholength=200/rholength = ${dir_shearloop_arr[$j]}/g" data_FR_Jing2.py 
    python3 modify_data.py  
    #sed -i "s/rholength=200/rholength = ${dir_shearloop_arr[$j]}/g" ${dir_vloop_arr[$i]}_${dir_shearloop_arr[$j]}_loop.ctrl
    ./paradisJHU ${dir_vloop_arr[$i]}_${dir_shearloop_arr[$j]}_loop.ctrl
    echo $i
#        filename = ${mod_fr_data_arr[$i]}'_'${dir_shearloop_arr[$j]}'_frank_read.data'
#        scp $filename ygu23@jhu.edu@login.marcc.jhu.edu:/home-4/ygu23@jhu.edu/work/yejun/AM/Composite/"$cell_name$counter/Prerun/frank_read.data"
    cd ../..
done
#for dir in '/home/yejun/Desktop/Works/Composite/'$dir_seq$'*/'     # list directories in the form 
#do
#    dir=${dir%*/}      # remove the trailing "/"
#    echo ${dir##*/}    # print everything after the final "/"
#done
