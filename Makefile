.PHONY: clean All

All:
	@echo "----------Building project:[ primesTo2Bil - Debug ]----------"
	@"$(MAKE)" -f  "primesTo2Bil.mk"
clean:
	@echo "----------Cleaning project:[ primesTo2Bil - Debug ]----------"
	@"$(MAKE)" -f  "primesTo2Bil.mk" clean
