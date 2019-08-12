##
## Auto Generated makefile by CodeLite IDE
## any manual changes will be erased      
##
## Debug
ProjectName            :=Assignment
ConfigurationName      :=Debug
WorkspacePath          :=/home/jl5916/Desktop/HPC
ProjectPath            :=/home/jl5916/Desktop/HPC/Assignment
IntermediateDirectory  :=./Debug
OutDir                 := $(IntermediateDirectory)
CurrentFileName        :=
CurrentFilePath        :=
CurrentFileFullPath    :=
User                   :=Javier Leguina Peral
Date                   :=14/03/19
CodeLitePath           :=/home/jl5916/.codelite
LinkerName             :=g++
SharedObjectLinkerName :=g++ -shared -fPIC
ObjectSuffix           :=.o
DependSuffix           :=.o.d
PreprocessSuffix       :=.o.i
DebugSwitch            :=-gstab
IncludeSwitch          :=-I
LibrarySwitch          :=-l
OutputSwitch           :=-o 
LibraryPathSwitch      :=-L
PreprocessorSwitch     :=-D
SourceSwitch           :=-c 
OutputFile             :=$(IntermediateDirectory)/$(ProjectName)
Preprocessors          :=
ObjectSwitch           :=-o 
ArchiveOutputSwitch    := 
PreprocessOnlySwitch   :=-E 
ObjectsFileList        :="Assignment.txt"
PCHCompileFlags        :=
MakeDirCommand         :=mkdir -p
LinkOptions            :=  
IncludePath            :=  $(IncludeSwitch). $(IncludeSwitch). 
IncludePCH             := 
RcIncludePath          := 
Libs                   := 
ArLibs                 :=  
LibPath                := $(LibraryPathSwitch). 

##
## Common variables
## AR, CXX, CC, AS, CXXFLAGS and CFLAGS can be overriden using an environment variables
##
AR       := ar rcus
CXX      := g++
CC       := gcc
CXXFLAGS :=  -g -O0 -Wall $(Preprocessors)
CFLAGS   :=  -g -O0 -Wall $(Preprocessors)
ASFLAGS  := 
AS       := as


##
## User defined environment variables
##
CodeLiteDir:=/usr/share/codelite
Objects0=$(IntermediateDirectory)/burguers.cpp$(ObjectSuffix) $(IntermediateDirectory)/Model.cpp$(ObjectSuffix) $(IntermediateDirectory)/Burgers.cpp$(ObjectSuffix) 



Objects=$(Objects0) 

##
## Main Build Targets 
##
.PHONY: all clean PreBuild PrePreBuild PostBuild MakeIntermediateDirs
all: $(OutputFile)

$(OutputFile): $(IntermediateDirectory)/.d $(Objects) 
	@$(MakeDirCommand) $(@D)
	@echo "" > $(IntermediateDirectory)/.d
	@echo $(Objects0)  > $(ObjectsFileList)
	$(LinkerName) $(OutputSwitch)$(OutputFile) @$(ObjectsFileList) $(LibPath) $(Libs) $(LinkOptions)

MakeIntermediateDirs:
	@test -d ./Debug || $(MakeDirCommand) ./Debug


$(IntermediateDirectory)/.d:
	@test -d ./Debug || $(MakeDirCommand) ./Debug

PreBuild:


##
## Objects
##
$(IntermediateDirectory)/burguers.cpp$(ObjectSuffix): burguers.cpp $(IntermediateDirectory)/burguers.cpp$(DependSuffix)
	$(CXX) $(IncludePCH) $(SourceSwitch) "/home/jl5916/Desktop/HPC/Assignment/burguers.cpp" $(CXXFLAGS) $(ObjectSwitch)$(IntermediateDirectory)/burguers.cpp$(ObjectSuffix) $(IncludePath)
$(IntermediateDirectory)/burguers.cpp$(DependSuffix): burguers.cpp
	@$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) -MG -MP -MT$(IntermediateDirectory)/burguers.cpp$(ObjectSuffix) -MF$(IntermediateDirectory)/burguers.cpp$(DependSuffix) -MM burguers.cpp

$(IntermediateDirectory)/burguers.cpp$(PreprocessSuffix): burguers.cpp
	$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) $(IntermediateDirectory)/burguers.cpp$(PreprocessSuffix) burguers.cpp

$(IntermediateDirectory)/Model.cpp$(ObjectSuffix): Model.cpp $(IntermediateDirectory)/Model.cpp$(DependSuffix)
	$(CXX) $(IncludePCH) $(SourceSwitch) "/home/jl5916/Desktop/HPC/Assignment/Model.cpp" $(CXXFLAGS) $(ObjectSwitch)$(IntermediateDirectory)/Model.cpp$(ObjectSuffix) $(IncludePath)
$(IntermediateDirectory)/Model.cpp$(DependSuffix): Model.cpp
	@$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) -MG -MP -MT$(IntermediateDirectory)/Model.cpp$(ObjectSuffix) -MF$(IntermediateDirectory)/Model.cpp$(DependSuffix) -MM Model.cpp

$(IntermediateDirectory)/Model.cpp$(PreprocessSuffix): Model.cpp
	$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) $(IntermediateDirectory)/Model.cpp$(PreprocessSuffix) Model.cpp

$(IntermediateDirectory)/Burgers.cpp$(ObjectSuffix): Burgers.cpp $(IntermediateDirectory)/Burgers.cpp$(DependSuffix)
	$(CXX) $(IncludePCH) $(SourceSwitch) "/home/jl5916/Desktop/HPC/Assignment/Burgers.cpp" $(CXXFLAGS) $(ObjectSwitch)$(IntermediateDirectory)/Burgers.cpp$(ObjectSuffix) $(IncludePath)
$(IntermediateDirectory)/Burgers.cpp$(DependSuffix): Burgers.cpp
	@$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) -MG -MP -MT$(IntermediateDirectory)/Burgers.cpp$(ObjectSuffix) -MF$(IntermediateDirectory)/Burgers.cpp$(DependSuffix) -MM Burgers.cpp

$(IntermediateDirectory)/Burgers.cpp$(PreprocessSuffix): Burgers.cpp
	$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) $(IntermediateDirectory)/Burgers.cpp$(PreprocessSuffix) Burgers.cpp


-include $(IntermediateDirectory)/*$(DependSuffix)
##
## Clean
##
clean:
	$(RM) -r ./Debug/


