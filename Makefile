SRC = options.cu determ.cu velocity.cu collision_integral.cu compute_matrix.cu velocity_grid.cu
OBJECTS=$(SRC:.cu=.o)
OBJDIR = obj
EXECUTABLE = determ
OBJLIST = $(addprefix $(OBJDIR)/,$(OBJECTS))

rebuild: 
	$(MAKE) clean 
	$(MAKE) $(OBJDIR)
	$(MAKE) $(EXECUTABLE)
#all: $(OBJDIR) $(EXECUTABLE) 

$(OBJDIR):
	mkdir $(OBJDIR)

$(OBJDIR)/%.o: %.cu
	nvcc -c $< -o $@

$(EXECUTABLE): | $(OBJLIST)
	nvcc $(OBJLIST) -o $@

clean: 
	rm -f obj/* $(EXECUTABLE)

-include $(DEPS)
