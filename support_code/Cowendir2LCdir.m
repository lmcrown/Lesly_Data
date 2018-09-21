function [LCpath]=Cowendir2LCdir(cowen_path)

[path,name,ext]=fileparts(cowen_path);

LCpath=[pwd '\' name ext];
