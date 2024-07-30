Input: Run Settings

Examples in DIV3D/test/…/run_settings.nml

Section: bfield_nml

- rmp_type = specifies format of magnetic field info
- bgrid_fname = bgrid path + file name (if rmp_type=bgrid)
- vmec_coils_file = path + file name of VMEC coils file
- vmec_extcur_set = external current scaling

Section: run_settings

- fname_plist = parts.list file name
- fname_ves = vessel.part file name + path
- fname_surf = surface line output file name
- fname_launch = launch points output file name
- fname_parts = all parts output file name
- fname_hit = fieldline trace prior to strikepoint output file name
- fname_intpts = intersection points output file name
- fname_nhit = number of hits output file name
- fname_ptri = part triangular groups output file name
- fname_ptri_mid = triangular mid-points output file name
- nfp = number of field periods
- Rstart, Zstart, Phistart = location of LCFS trace start point
- dphi_line_surf_deg = degree of accuracy, per integration step for LCFS surface field line tracing without diffusion
- ntran_surf = number of toroidal transits (make smaller to make surf_line.out file size smaller)
- npts_start = number of points randomly distributed along the surface line
- dmag = diffusion coefficient in m2 / m
- dphi_line_diff_deg = degree of accuracy, per integration step for field line tracing with diffusion
- ntran_diff = max number of toroidal transits
- trace_surface_opt = whether to trace out surface, false means read in existing surface file
- myseed = ? random number generator seed?
- hit_length = length of end of fieldline recorded
- lsfi_tol = tolerance of adjacent facets aligning

---

Input: Magnetic Field Grid

Magnetic field options:

- vmec_coils – self explanatory. This is probably relatively slow
- vmec_coils_to_fil – This probably turns the VMEC coilset into a series of current filaments. Probably I was just trying for the most efficient method here.
- Xdr – Format used by W7X, basically just RZ cartesian slices with the RZ grid shifting with the axis toroidally
- Bgrid – Just a RZ cartesian grid at some toroidal resolution. Most general and can be very fast.

---

Input: Parts List

Examples in DIV3D/test/…/parts.list

Lists out .part files and their path locations

