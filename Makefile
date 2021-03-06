############################################
#############################################
##          Makefile for Knoblauch
#############################################
#############################################

CXX = g++

#############################################
##     Name of executable
#############################################
EXECUTABLE = Knoblauch

#############################################
##     Name of "tester" executable
#############################################
TESTER = Schnittlauch

#############################################
##   Dicrectories
#############################################
DEPEND_DIR = .DEPEND/

VERSION_DIR = Version/

INPUT_DIR = Input/

LINALG_DIR = Lin_Alg/

PHYSICS_DIR = Physics/

GRID1D_DIR = Grid1D/

GRID2D_DIR = Grid2D/

DISCRETIZATIONS2D_DIR = Discretizations2D/

############################################
############################################
##  Compilation related stuff
############################################
############################################
EIGEN_INCLUDE = /usr/include/eigen3
VTK_INCLUDE   = /usr/include/vtk
VTK_LIBRARY   = /usr/lib64/vtk

CXXFLAGS = -O3 -DNDEBUG -Wall -Wextra -Wpedantic -std=c++11 \
           -isystem $(EIGEN_INCLUDE) -isystem $(VTK_INCLUDE) -Wno-deprecated

DEPEND_FLAGS = $(CXXFLAGS) -MM -MF $(addprefix $(DEPEND_DIR),$*.d) -MT $*.o -MP

LIBS = -L$(VTK_LIBRARY) -lvtkCommonCore -lvtkCommonDataModel  -lvtkIOCore -lvtkIOXML
#LIBS = -L$(VTK_LIBRARY) -lvtkCommon -lvtkIO -lvtkFiltering
TESTER_LIBS = -lgtest -lm

############################################
############################################
##     "Special" Version-related stuff
############################################
############################################
PYTHON = $(shell which python3)

VERSION_SCRIPT = Version_Info_Generator.py

GIT_VERSION := $(shell $(PYTHON) $(VERSION_DIR)$(VERSION_SCRIPT))

############################################
############################################
##          Source files
############################################
############################################
SRC_KNOBLAUCH = Knoblauch_main.cpp
SRC_SCHNITTLAUCH = Schnittlauch_main.cpp

############################################
##   Version info
SRC_VERSION = Name_Info.cpp $(GIT_VERSION)
SRC_VERSION_TESTS = Git_Version_Tests.cpp

SRC_VERSION := $(addprefix $(VERSION_DIR), $(SRC_VERSION))
SRC_VERSION_TESTS := $(addprefix $(VERSION_DIR), $(SRC_VERSION_TESTS))

############################################
##  Input-Parameter Handling
SRC_INPUT = Input_Parameter_Tree.cpp
SRC_INPUT_TESTS = Input_Helpers_Tests.cpp

SRC_INPUT := $(addprefix $(INPUT_DIR), $(SRC_INPUT))
SRC_INPUT_TESTS := $(addprefix $(INPUT_DIR), $(SRC_INPUT_TESTS))

############################################
##   Linear Algebra
SRC_LINALG =
SRC_LINALG_TESTS = Basic_Linear_Algebra_Tests.cpp

SRC_LINALG := $(addprefix $(LINALG_DIR), $(SRC_LINALG))
SRC_LINALG_TESTS := $(addprefix $(LINALG_DIR), $(SRC_LINALG_TESTS))

############################################
##   Physics
SRC_PHYSICS =
SRC_PHYSICS_TESTS = Calorically_Perfect_Gas_Tests.cpp

SRC_PHYSICS := $(addprefix $(PHYSICS_DIR), $(SRC_PHYSICS))
SRC_PHYSICS_TESTS := $(addprefix $(PHYSICS_DIR), $(SRC_PHYSICS_TESTS))

############################################
##   Grid1D
SRC_GRID1D =
SRC_GRID1D_TESTS = Block_Uniform_1D_Tests.cpp

SRC_GRID1D := $(addprefix $(GRID1D_DIR), $(SRC_GRID1D))
SRC_GRID1D_TESTS := $(addprefix $(GRID1D_DIR), $(SRC_GRID1D_TESTS))

############################################
##   Grid2D
SRC_GRID2D = Make_Cartesian_Block_2D.cpp \
             Make_Cylinder_Block_2D.cpp \
             Write_Output_VTK.cpp \
             Grid_Tree/Refine_Block2D_I.cpp \
             Grid_Tree/Refine_Block2D_J.cpp \
             Grid_Tree/Refine_Block2D_IJ.cpp

SRC_GRID2D_TESTS = Geometry_2D_Tests.cpp \
                   Block_2D_Tests.cpp \
                   Cartesian_Block_2D_Tests.cpp \
                   Make_Cylinder_Block_2D_Tests.cpp \
                   Curve2D/Line_Segment2D_Tests.cpp \
                   Curve2D/Circular_Arc2D_Tests.cpp \
                   Curve2D/Spline2D_Tests.cpp

SRC_GRID2D := $(addprefix $(GRID2D_DIR), $(SRC_GRID2D))
SRC_GRID2D_TESTS := $(addprefix $(GRID2D_DIR), $(SRC_GRID2D_TESTS))

############################################
##   SpatialDiscretizations2D
SRC_DISCRETIZATIONS2D =
SRC_DISCRETIZATIONS2D_TESTS =

SRC_DISCRETIZATIONS2D := \
      $(addprefix $(DISCRETIZATIONS2D_DIR), $(SRC_DISCRETIZATIONS2D))
SRC_DISCRETIZATIONS2D_TESTS := \
      $(addprefix $(DISCRETIZATIONS2D_DIR), $(SRC_DISCRETIZATIONS2D_TESTS))

############################################
##   All sources for Knbolauch
SRC_MAIN = $(SRC_VERSION) \
           $(SRC_INPUT) \
           $(SRC_LINALG) \
           $(SRC_PHYSICS) \
           $(SRC_GRID1D) \
           $(SRC_GRID2D) \
           $(SRC_DISCRETIZATIONS2D)

############################################
##   All sources for Schnittlauch
SRC_TESTER = $(SRC_MAIN) \
             $(SRC_VERSION_TESTS) \
             $(SRC_INPUT_TESTS) \
             $(SRC_LINALG_TESTS) \
             $(SRC_PHYSICS_TESTS) \
             $(SRC_GRID1D_TESTS) \
             $(SRC_GRID2D_TESTS) \
             $(SRC_DISCRETIZATIONS2D_TESTS)

############################################
##   All sources
SRC_ALL = $(SRC_TESTER) \
          $(SRC_KNOBLAUCH) \
          $(SRC_SCHNITTLAUCH)

############################################
############################################
##          Object files
############################################
############################################
OBJ_KNOBLAUCH    = $(SRC_KNOBLAUCH:.cpp=.o)
OBJ_SCHNITTLAUCH = $(SRC_SCHNITTLAUCH:.cpp=.o)

OBJ_MAIN         = $(SRC_MAIN:.cpp=.o)
OBJ_TESTER       = $(SRC_TESTER:.cpp=.o)

OBJ_ALL          = $(SRC_ALL:.cpp=.o)

############################################
############################################
##           Targets
############################################
############################################
ifndef MAKECMDGOALS
MAKECMDGOALS = $(EXECUTABLE)
endif

PHONY_TARGETS = all clean new out test

.phony : $(PHONY_TARGETS)

$(EXECUTABLE) : $(OBJ_MAIN) $(OBJ_KNOBLAUCH)
	$(CXX) $(CXXFLAGS) $(OBJ_MAIN) $(OBJ_KNOBLAUCH) -o $(EXECUTABLE) $(LIBS)

$(TESTER) : $(OBJ_TESTER) $(OBJ_SCHNITTLAUCH)
	$(CXX) $(CXXFLAGS) $(OBJ_TESTER) $(OBJ_SCHNITTLAUCH) -o $(TESTER) $(LIBS) $(TESTER_LIBS)

test : $(TESTER)
	@./$(TESTER)

all : $(EXECUTABLE) $(TESTER)

clean :
	rm -f $(OBJ_ALL)
	rm -fr $(DEPEND_DIR)

new : clean
	rm -f gmon.out
	rm -f core*
	rm -f $(EXECUTABLE)
	rm -f $(TESTER)

out :
	@echo "Smooch!"

############################################
############################################
##          Dependencies
############################################
############################################
$(DEPEND_DIR)%.d : %.cpp
	@echo "Generating dependency information for$(patsubst %.cpp, %.o, $<)"
	@mkdir -p $(dir $@)
	@$(CXX) -c $(DEPEND_FLAGS) $<
	@sed -i 's@^\(.*\.o\):@$@ \1: Makefile@' $@

NO_DEP_TARGETS = clean new out

ifneq ($(strip $(filter-out $(NO_DEP_TARGETS), $(MAKECMDGOALS)) ), )
-include $(patsubst %.cpp,$(DEPEND_DIR)%.d, $(SRC_ALL))
endif
