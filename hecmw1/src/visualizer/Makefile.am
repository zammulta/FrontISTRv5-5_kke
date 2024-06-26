MPIDIR            = @mpidir@
MPIBINDIR         = @mpibindir@
MPILIBDIR         = @mpilibdir@
MPIINCDIR         = @mpiincdir@
MPILIBS           = @mpilibs@

PREFIX            = @prefix@
BINDIR            = @bindir@
LIBDIR            = @libdir@
INCLUDEDIR        = @includedir@

HECMWBINDIR       = ../../bin
HECMWLIBDIR       = ../../lib
HECMWINCDIR       = ../../include
HECMWLIBS         = @hecmwlibs@

REFINERDIR        = @refinerdir@
REFINERINCDIR     = @refinerincdir@
REFINERLIBDIR     = @refinerlibdir@
REFINERLIBS       = @refinerlibs@

CC                = @cc@
CFLAGS            = @cflags@
BASE_CFLAGS       = @base_cflags@
MPI_CFLAGS        = @mpi_cflags@
HECMW_CFLAGS      = @hecmw_cflags@
ALL_CFLAGS        = $(BASE_CFLAGS) $(HECMW_CFLAGS) $(MPI_CFLAGS) $(CFLAGS)
LDFLAGS           = @ldflags@
MPI_LDFLAGS       = @mpi_ldflags@
HECMW_LDFLAGS     = @hecmw_ldflags@
ALL_LDFLAGS       = $(HECMW_LDFLAGS) $(MPI_LDFLAGS) $(LDFLAGS)
OPTFLAGS          = @optflags@

F90               = @f90@
F90FLAGS          = @f90flags@
BASE_F90FLAGS     = @base_f90flags@
MPI_F90FLAGS      = @mpi_f90flags@
HECMW_F90FLAGS    = @hecmw_f90flags@
ALL_F90FLAGS      = $(BASE_F90FLAGS) $(HECMW_F90FLAGS) $(MPI_F90FLAGS) $(F90FLAGS)
F90LDFLAGS        = @f90ldflags@
MPI_F90LDFLAGS    = @mpi_f90ldflags@
HECMW_F90LDFLAGS  = @hecmw_f90ldflags@
ALL_F90LDFLAGS    = $(HECMW_F90LDFLAGS) $(MPI_F90LDFLAGS) $(F90LDFLAGS)
F90OPTFLAGS       = @f90optflags@

AR                = @ar@
CP                = @cp@
RM                = @rm@
RANLIB            = @ranlib@
MKDIR             = @mkdir@

#-----------------------------------------------------------------------------#
TARGET            = @hecmwlib_targetfile@
F90TARGET         = @hecmwlib_f90targetfile@

F90MODULEPOSTFIX  = @f90modulepostfix@
COBJFILEPOSTFIX   = @cobjfilepostfix@
F90OBJFILEPOSTFIX = @f90objfilepostfix@

#-----------------------------------------------------------------------------#
.SUFFIXES:
.SUFFIXES: .@cobjfilepostfix@ .c

.c.@cobjfilepostfix@:
	$(CC) -c $< $(ALL_CFLAGS) $(OPTFLAGS)

#-----------------------------------------------------------------------------#
PSF_OBJS = \
	hecmw_vis_surface_main.@cobjfilepostfix@ \
	hecmw_vis_comm_util.@cobjfilepostfix@ \
	hecmw_vis_mem_util.@cobjfilepostfix@ \
	hecmw_vis_connectivity_build.@cobjfilepostfix@ \
	hecmw_vis_surface_compute.@cobjfilepostfix@ \
	hecmw_vis_case_table.@cobjfilepostfix@ \
	hecmw_vis_calc_attr.@cobjfilepostfix@ \
	hecmw_vis_intersection_find.@cobjfilepostfix@ \
	hecmw_vis_patch_const.@cobjfilepostfix@ \
	hecmw_vis_tetra_intersect.@cobjfilepostfix@ \
	hecmw_vis_combine.@cobjfilepostfix@ \
	hecmw_vis_rendering.@cobjfilepostfix@ \
	hecmw_vis_define_parameters.@cobjfilepostfix@ \
	hecmw_vis_font_texture.@cobjfilepostfix@ \
	hecmw_vis_color_mapping.@cobjfilepostfix@ \
	hecmw_vis_color_composite_sf.@cobjfilepostfix@ \
	hecmw_vis_generate_histogram_sf.@cobjfilepostfix@ \
	hecmw_vis_subimage_composite_sf.@cobjfilepostfix@ \
	hecmw_fstr_output_femap.@cobjfilepostfix@ \
	hecmw_fstr_output_vtk.@cobjfilepostfix@ \
	hecmw_fstr_endian.@cobjfilepostfix@

PVR_OBJS = \
	hecmw_vis_pvr_main.@cobjfilepostfix@ \
	hecmw_vis_comm_util.@cobjfilepostfix@ \
	hecmw_vis_mem_util.@cobjfilepostfix@ \
	hecmw_vis_new_refine.@cobjfilepostfix@ \
	hecmw_vis_ray_trace.@cobjfilepostfix@ \
	hecmw_vis_define_parameters.@cobjfilepostfix@ \
	hecmw_vis_font_texture.@cobjfilepostfix@ \
	hecmw_vis_color_mapping.@cobjfilepostfix@ \
	hecmw_vis_ucd_trans.@cobjfilepostfix@ \
	hecmw_vis_voxel_gen.@cobjfilepostfix@ \
	hecmw_vis_read_voxel.@cobjfilepostfix@ \
	hecmw_vis_mesh_extent.@cobjfilepostfix@ \
	hecmw_vis_color_composite_vr.@cobjfilepostfix@ \
	hecmw_vis_generate_histogram_vr.@cobjfilepostfix@ \
	hecmw_vis_subimage_composite_vr.@cobjfilepostfix@

OBJS = \
	hecmw_vis_read_control.@cobjfilepostfix@ \
	hecmw_visualizer.@cobjfilepostfix@ \
	$(PSF_OBJS) \
	$(PVR_OBJS) \

HEADERS = \
	hecmw_fstr_output_femap.h \
	hecmw_vis_SF_geom.h \
	hecmw_vis_bmp.h \
	hecmw_vis_calc_attr.h \
	hecmw_vis_case_table.h \
	hecmw_vis_color_composite_sf.h \
	hecmw_vis_color_composite_vr.h \
	hecmw_vis_color_mapping.h \
	hecmw_vis_combine.h \
	hecmw_vis_comm_util.h \
	hecmw_vis_connectivity_build.h \
	hecmw_vis_define_parameters.h \
	hecmw_vis_define_parameters_vr.h \
	hecmw_vis_endian.h \
	hecmw_vis_font_texture.h \
	hecmw_vis_generate_histogram_sf.h \
	hecmw_vis_generate_histogram_vr.h \
	hecmw_vis_intersection_find.h \
	hecmw_vis_mem_util.h \
	hecmw_vis_mesh_extent.h \
	hecmw_vis_new_refine.h \
	hecmw_vis_patch_const.h \
	hecmw_vis_psf_rendering.h \
	hecmw_vis_pvr_main.h \
	hecmw_vis_ray_trace.h \
	hecmw_vis_read_control.h \
	hecmw_vis_read_voxel.h \
	hecmw_vis_rendering.h \
	hecmw_vis_resampling.h \
	hecmw_vis_subimage_composite_sf.h \
	hecmw_vis_subimage_composite_vr.h \
	hecmw_vis_surface_compute.h \
	hecmw_vis_surface_main.h \
	hecmw_vis_tetra_intersect.h \
	hecmw_vis_ucd_trans.h \
	hecmw_vis_voxel_gen.h \
	hecmw_visualizer.h

INSTALLHEADERS = \
        hecmw_vis_mem_util.h \
        hecmw_vis_comm_util.h \
        hecmw_vis_read_control.h \
        hecmw_vis_SF_geom.h \
        hecmw_vis_ray_trace.h \
        hecmw_vis_psf_rendering.h \
        hecmw_visualizer.h

#-----------------------------------------------------------------------------#
all: @build_target@

build-default: $(TARGET)

build-serial: $(TARGET)

build-without-f: $(TARGET)

$(TARGET): $(OBJS)
	$(AR) $(HECMWLIBDIR)/$@ $(OBJS)

$(PSF_OBJS): $(HEADERS)

$(PVR_OBJS): $(HEADERS)

install-header: $(HEADERS)
	$(CP) $(INSTALLHEADERS) $(HECMWINCDIR)/.

install-module:
	$(CP) *.$(F90MODULEPOSTFIX) $(HECMWINCDIR)/.

install:

clean:
	$(RM) *.$(COBJFILEPOSTFIX) *.$(F90OBJFILEPOSTFIX) *.$(F90MODULEPOSTFIX)

distclean: clean

depend:
