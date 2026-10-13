# ares_domain_overview.pml -- domain overview of aRES (UC5) with the Zn(Cys)4 ribbon.
#
# Usage:   pymol ares_domain_overview.pml
#   or inside PyMOL:   @analysis/scripts/ares_domain_overview.pml
#
# Set the path to your model below, then run. The script builds four scenes:
#   F1  structural (PAE) domains      <- the one to use
#   F2  Pfam annotation, for contrast
#   F3  coloured by pLDDT confidence
#   F4  zinc site close-up
#
# Boundaries come from docs/05 section 1.1 and docs/07 section 1. They are the
# AF3 PAE structural blocks, NOT the Pfam call -- the two disagree on purpose,
# see the F2 scene and the note at the bottom of this file.

delete all

# >>> EDIT THIS LINE: point it at your model <<<
load aRES_Zn_segA.pdb, aRES

# ---------------------------------------------------------------- preferences
bg_color white
set ray_opaque_background, 0
hide everything
set cartoon_transparency, 0
set cartoon_side_chain_helper, 1
set sphere_scale, 0.45
set dash_gap, 0.25
set dash_width, 2.5
set label_size, 16
set label_color, black
set antialias, 2
set ray_shadows, 0
set specular, 0.2
set valence, 0

# ------------------------------------------------- structural (PAE) selections
# NOTE: the four zinc-binding cysteines are typed CYF in this file, not CYS.
# A selection written as "resn CYS" silently misses all four.
select zn_ribbon,  aRES and resi 1-48
select ntd1,       aRES and resi 49-120
select hinge,      aRES and resi 121-125
select res_block,  aRES and resi 126-447
select flexloop,   aRES and resi 392-408
select zn_cys,     aRES and resn CYF
select zn_metal,   aRES and elem ZN

# Pfam calls, for the contrast scene
select pf18870,    aRES and resi 49-136
select pf08808,    aRES and resi 199-363

# landmarks that verify against the paper (docs/08 line 22)
select catalytic,  aRES and resi 182+184+215+236
# D104: putative gp5 Ser122 contact -- single AF3 cofold at poor ipTM, see docs/05 section 2
select d104,       aRES and resi 104
# free cysteines with no metal: oxidation risk, docs/08 section 3d
select free_cys,   aRES and resi 78+96+249

# NOT selected on purpose: residues 234 and 293. This model has Ile234/Ser293
# where the project's defence-null list reads S234A/V293A -- unresolved frame
# discrepancy, docs/08 section 3e. Do not label them until that is settled.

# =========================================================== F1 PAE domains
show cartoon, aRES
color grey70,    aRES
color deepteal,  zn_ribbon
color orange,    ntd1
color red,       hinge
color skyblue,   res_block
color palegreen, flexloop

show spheres, zn_metal
color grey40,  zn_metal
show sticks,   zn_cys and not name N+C+O
color yellow,  zn_cys and elem S
show sticks,   d104 and not name N+C+O
color magenta, d104 and elem O

orient aRES
scene F1_domains, store, \
  "aRES PAE domains: Zn ribbon 1-48 teal | NTD-1 49-120 orange | hinge 121-125 red | RES 126-447 blue"

# =========================================================== F2 Pfam contrast
disable all
enable aRES
color grey80, aRES
color orange, pf18870
color skyblue, pf08808
color red,    hinge
show spheres, zn_metal

# The point of this scene: PF18870 (49-136) runs straight through the hinge,
# and PF08808 (199-363) excludes R182/R184, which are defence-null. A construct
# built on the Pfam RES boundary would lack two active-site arginines.
show spheres, catalytic and name CA
color firebrick, catalytic
scene F2_pfam, store, \
  "Pfam: PF18870 49-136 crosses the hinge; PF08808 199-363 excludes R182/R184 (red spheres)"

# =========================================================== F3 pLDDT
# B-factor column in this file is a pLDDT-like confidence field: HIGH = confident.
# Verified against the per-region means recorded in docs/07.
disable all
enable aRES
hide everything, aRES
show cartoon, aRES
spectrum b, red_yellow_green, aRES, 40, 95
show spheres, zn_metal
color grey40, zn_metal
scene F3_plddt, store, \
  "pLDDT 40-95 (red=low). Hinge 121-125 and loop 392-408 are the low-confidence regions"

# =========================================================== F4 zinc site
disable all
enable aRES
hide everything, aRES
show cartoon, aRES
color grey90, aRES
set cartoon_transparency, 0.75, aRES

show sticks, zn_cys and not name N+C+O
color yellow, zn_cys and elem S
color grey50, zn_cys and elem C
show spheres, zn_metal
color grey30, zn_metal

# the four coordination bonds -- measured at 2.242-2.252 A, mean 2.248
distance zn_bonds, zn_metal, zn_cys and name SG, 3.0
color black, zn_bonds
hide labels, zn_bonds

label zn_cys and name CA, "%s%s" % (resn, resi)
orient zn_cys
zoom zn_cys, 4
scene F4_zincsite, store, \
  "Zn(Cys)4: C8/C11/C31/C34, Zn-S 2.248 A mean, S-Zn-S RMS 3.87 deg from tetrahedral"

# ----------------------------------------------------------------------------
scene F1_domains, recall

print ""
print "  scenes: F1_domains  F2_pfam  F3_plddt  F4_zincsite"
print "  ray 1600, 1200   then   png overview.png, dpi=300"
print ""
print "  NOTE: this model is residues 1-447. The disordered tail 448-459 is"
print "        absent, so it cannot be shown. Full length is 1-459."
print ""
