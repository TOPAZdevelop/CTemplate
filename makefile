Here = $(PWD)
PDFDir = $(Here)/PDFs
ObjectDir = $(Here)/objects


CC = icc
FC = ifort

Exec = ./Integrate



VegasObj = $(ObjectDir)/nvegas.o
PSObj = $(ObjectDir)/genps.o
PDFObj = $(ObjectDir)/alphaS.o \
         $(ObjectDir)/mstwpdf.o
MainObj = $(ObjectDir)/integrate_crosssection.o




all:  $(general) $(VegasObj) $(PSObj) $(PDFObj) $(MainObj) 
	@echo " linking"
	@echo " executable file is " $(Exec)
	@echo " "
	$(CC) $(MainObj) $(PDFObj) $(PSObj) $(VegasObj) -lm -lifcore -o $(Exec)

clean:
	rm -f ./objects/*.o
	rm -f ./$(Exec)



general: *.c
	@echo " creating directories"
	mkdir -p $(Here)/objects 
	mkdir -p $(Here)/data


$(VegasObj): $(Here)/nvegas.c $(Here)/vegas.h
	@echo " compiling vegas" 
	$(CC) -c $< -o $@

$(PSObj): $(Here)/genps.c $(Here)/genps.h
	@echo " compiling phase space generator" 
	$(CC) -c $< -o $@

$(ObjectDir)/alphaS.o: $(PDFDir)/alphaS.f
	@echo " compiling alphaS" 
	$(FC) -c $< -o $@
	
$(ObjectDir)/mstwpdf.o: $(PDFDir)/mstwpdf.f
	@echo " compiling mstwpdf" 
	$(FC) -c $< -o $@

$(MainObj): $(Here)/integrate_crosssection.c
	@echo " compiling main objects" 
	$(CC) -c $< -o $@





# supresses command calls
.SILENT:
