##
## Auto Generated makefile by CodeLite IDE
## any manual changes will be erased      
##
## Debug
ProjectName            :=Chapter21
ConfigurationName      :=Debug
WorkspacePath          :="/home/jl5916/Desktop/High Performance Computing"
ProjectPath            :="/home/jl5916/Desktop/High Performance Computing/Chapter21"
IntermediateDirectory  :=./Debug
OutDir                 := $(IntermediateDirectory)
CurrentFileName        :=
CurrentFilePath        :=
CurrentFileFullPath    :=
User                   :=Javier Leguina Peral
Date                   :=04/03/19
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
ObjectsFileList        :="Chapter21.txt"
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
Objects0=$(IntermediateDirectory)/Ex21_2_generate.cpp$(ObjectSuffix) $(IntermediateDirectory)/Ex21_2_fftw.cpp$(ObjectSuffix) 



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
$(IntermediateDirectory)/Ex21_2_generate.cpp$(ObjectSuffix): Ex21_2/generate.cpp $(IntermediateDirectory)/Ex21_2_generate.cpp$(DependSuffix)
	$(CXX) $(IncludePCH) $(SourceSwitch) "/home/jl5916/Desktop/High Performance Computing/Chapter21/Ex21_2/generate.cpp" $(CXXFLAGS) $(ObjectSwitch)$(IntermediateDirectory)/Ex21_2_generate.cpp$(ObjectSuffix) $(IncludePath)
$(IntermediateDirectory)/Ex21_2_generate.cpp$(DependSuffix): Ex21_2/generate.cpp
	@$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) -MG -MP -MT$(IntermediateDirectory)/Ex21_2_generate.cpp$(ObjectSuffix) -MF$(IntermediateDirectory)/Ex21_2_generate.cpp$(DependSuffix) -MM Ex21_2/generate.cpp

$(IntermediateDirectory)/Ex21_2_generate.cpp$(PreprocessSuffix): Ex21_2/generate.cpp
	$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) $(IntermediateDirectory)/Ex21_2_generate.cpp$(PreprocessSuffix) Ex21_2/generate.cpp

$(IntermediateDirectory)/Ex21_2_fftw.cpp$(ObjectSuffix): Ex21_2/fftw.cpp $(IntermediateDirectory)/Ex21_2_fftw.cpp$(DependSuffix)
	$(CXX) $(IncludePCH) $(SourceSwitch) "/home/jl5916/Desktop/High Performance Computing/Chapter21/Ex21_2/fftw.cpp" $(CXXFLAGS) $(ObjectSwitch)$(IntermediateDirectory)/Ex21_2_fftw.cpp$(ObjectSuffix) $(IncludePath)
$(IntermediateDirectory)/Ex21_2_fftw.cpp$(DependSuffix): Ex21_2/fftw.cpp
	@$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) -MG -MP -MT$(IntermediateDirectory)/Ex21_2_fftw.cpp$(ObjectSuffix) -MF$(IntermediateDirectory)/Ex21_2_fftw.cpp$(DependSuffix) -MM Ex21_2/fftw.cpp

$(IntermediateDirectory)/Ex21_2_fftw.cpp$(PreprocessSuffix): Ex21_2/fftw.cpp
	$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) $(IntermediateDirectory)/Ex21_2_fftw.cpp$(PreprocessSuffix) Ex21_2/fftw.cpp


-include $(IntermediateDirectory)/*$(DependSuffix)
##
## Clean
##
clean:
	$(RM) -r ./Debug/


