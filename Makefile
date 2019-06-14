.PHONY: clean All

All:
	@echo "----------Building project:[ primesTo2Bil - Release ]----------"
	@"$(MAKE)" -f  "primesTo2Bil.mk"
clean:
	@echo "----------Cleaning project:[ primesTo2Bil - Release ]----------"
	@"$(MAKE)" -f  "primesTo2Bil.mk" clean
