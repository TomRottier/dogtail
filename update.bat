@echo off

@REM Delete .in and .dir files
del .\model\dogtail.dir
del .\model\dogtail.in

@REM Run Julia script to update equations
julia getequations.jl