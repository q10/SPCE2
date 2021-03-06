#
# Generated Makefile - do not edit!
#
# Edit the Makefile in the project folder instead (../Makefile). Each target
# has a -pre and a -post target defined where you can add customized code.
#
# This makefile implements configuration specific macros and targets.


# Environment
MKDIR=mkdir
CP=cp
GREP=grep
NM=nm
CCADMIN=CCadmin
RANLIB=ranlib
CC=gcc
CCC=g++
CXX=g++
FC=gfortran
AS=as

# Macros
CND_PLATFORM=MinGW-Windows
CND_CONF=Release
CND_DISTDIR=dist
CND_BUILDDIR=build

# Include project Makefile
include Makefile2

# Object Directory
OBJECTDIR=${CND_BUILDDIR}/${CND_CONF}/${CND_PLATFORM}

# Object Files
OBJECTFILES= \
	${OBJECTDIR}/src/energy.o \
	${OBJECTDIR}/src/umbrella_spce_hamiltonian.o \
	${OBJECTDIR}/src/main.o \
	${OBJECTDIR}/src/random.o \
	${OBJECTDIR}/src/water_system.o \
	${OBJECTDIR}/src/runtime.o \
	${OBJECTDIR}/src/ipair_distance_sampler.o \
	${OBJECTDIR}/src/rdf_sampler.o \
	${OBJECTDIR}/src/spce_hamiltonian.o \
	${OBJECTDIR}/src/rotation.o \
	${OBJECTDIR}/src/config_reader.o \
	${OBJECTDIR}/src/water.o \
	${OBJECTDIR}/src/ion.o \
	${OBJECTDIR}/src/ewald.o


# C Compiler Flags
CFLAGS=

# CC Compiler Flags
CCFLAGS=
CXXFLAGS=

# Fortran Compiler Flags
FFLAGS=

# Assembler Flags
ASFLAGS=

# Link Libraries and Options
LDLIBSOPTIONS=

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS}
	"${MAKE}"  -f nbproject/Makefile-${CND_CONF}.mk ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/spce2.exe

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/spce2.exe: ${OBJECTFILES}
	${MKDIR} -p ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}
	${LINK.cc} -o ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/spce2 ${OBJECTFILES} ${LDLIBSOPTIONS} 

${OBJECTDIR}/src/energy.o: src/energy.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/energy.o src/energy.cpp

${OBJECTDIR}/src/umbrella_spce_hamiltonian.o: src/umbrella_spce_hamiltonian.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/umbrella_spce_hamiltonian.o src/umbrella_spce_hamiltonian.cpp

${OBJECTDIR}/src/main.o: src/main.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/main.o src/main.cpp

${OBJECTDIR}/src/random.o: src/random.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/random.o src/random.cpp

${OBJECTDIR}/src/water_system.o: src/water_system.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/water_system.o src/water_system.cpp

${OBJECTDIR}/src/runtime.o: src/runtime.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/runtime.o src/runtime.cpp

${OBJECTDIR}/src/ipair_distance_sampler.o: src/ipair_distance_sampler.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/ipair_distance_sampler.o src/ipair_distance_sampler.cpp

${OBJECTDIR}/src/rdf_sampler.o: src/rdf_sampler.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/rdf_sampler.o src/rdf_sampler.cpp

${OBJECTDIR}/src/spce_hamiltonian.o: src/spce_hamiltonian.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/spce_hamiltonian.o src/spce_hamiltonian.cpp

${OBJECTDIR}/src/rotation.o: src/rotation.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/rotation.o src/rotation.cpp

${OBJECTDIR}/src/config_reader.o: src/config_reader.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/config_reader.o src/config_reader.cpp

${OBJECTDIR}/src/water.o: src/water.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/water.o src/water.cpp

${OBJECTDIR}/src/ion.o: src/ion.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/ion.o src/ion.cpp

${OBJECTDIR}/src/ewald.o: src/ewald.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/ewald.o src/ewald.cpp

# Subprojects
.build-subprojects:

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r ${CND_BUILDDIR}/${CND_CONF}
	${RM} ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/spce2.exe

# Subprojects
.clean-subprojects:

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
