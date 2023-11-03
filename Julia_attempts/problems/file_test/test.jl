cd(@__DIR__)
pwd()
touch("but.txt")

f = open("but.txt","r")
close(f)
open("but.txt") do file
    for li in eachline(file) 
        println("$(li)") 
    end 
end 



using DelimitedFiles 
using Dates
using Plots

flname=Dates.format(now(), "YYYY_mm_dd_HH_MM") *"_test.csv"
touch(flname)
A =[[0:5...] [5:10...]]
#Writing contents to the file 
open(flname, "w") do io 
    writedlm( flname,  A, ',')
end; 
  
#readlm() method to read the DelimitedFiles 
tst = readdlm(flname,',')

plot(tst)

## ok now two DelimitedFiles

fld = Dates.format(now(), "YYYY_mm_dd_HH_MM")
mkdir(fld)
flname="test_1.csv"
touch(fld*"/"*flname)

A =[[0:5...] [5:10...]]
open(fld*"/"*flname, "w") do io 
    writedlm( fld*"/"*flname,  A, ',')
end; 
tst = readdlm(fld*"/"*flname,',')

plot(tst)
flname="test_2.csv"
touch(fld*"/"*flname)

B =[[10:15...] [15:20...]]
#Writing contents to the file 
open(fld*"/"*flname, "w") do io 
    writedlm( fld*"/"*flname,  B, ',')
end; 
  
#readlm() method to read the DelimitedFiles 
tst = readdlm(fld*"/"*flname,',')

plot(tst)

files = readdir(fld)

arr1 = []
arr2 = []
for fl in files
    tst = transpose(readdlm(fld*"/"*fl,','))
    append!(arr1,tst[1,:])
    append!(arr2,tst[2,:])
end

