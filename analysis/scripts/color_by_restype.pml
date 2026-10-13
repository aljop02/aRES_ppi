# color_by_restype.pml -- model-free charge map of a structure surface.
#
# Colours the surface by residue identity only: Lys/Arg blue, Asp/Glu red,
# His pale blue, everything else grey. No charge model, no protonation
# assumptions, no dependence on how the Zn site was parameterised -- which is
# the point. Use it as the control on an APBS / vacuum_electrostatics map.
#
# Usage:  pymol model.pdb color_by_restype.pml

bg_color white
hide everything
show surface, polymer
color grey70,    polymer
color red,       polymer and resn ASP+GLU
color blue,      polymer and resn LYS+ARG
color skyblue,   polymer and resn HIS
show spheres,    elem ZN
color magenta,   elem ZN

set surface_quality, 1
set ray_shadows, 0
set antialias, 2
set transparency, 0

# The two N-block modules, for reference.
select ribbon48, polymer and resi 1-48
select ntd1,     polymer and resi 49-120
deselect

# Front and back of the N-block.
orient polymer and resi 1-120
ray 1600, 1200
png figs/restype_front.png, dpi=300
turn y, 180
ray 1600, 1200
png figs/restype_back.png, dpi=300
turn y, 180

# Artifact control: same electrostatics calculation without the metal.
# If the KRKRK face is still blue on "nozn", the basic patch is real and the
# zinc was only adding gloss.
create nozn, polymer
