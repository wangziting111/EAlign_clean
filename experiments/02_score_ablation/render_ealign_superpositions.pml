set quiet, 1
set suspend_updates, 1
bg_color white
set ray_opaque_background, 0
set ray_shadows, 0
set orthoscopic, 1
set antialias, 2
set ray_trace_mode, 1
set ambient, 0.46
set direct, 0.58
set specular, 0.28
set shininess, 24
set depth_cue, 0
set cartoon_fancy_helices, 1
set cartoon_smooth_loops, 1
set cartoon_flat_sheets, 1
set cartoon_sampling, 14
set cartoon_loop_radius, 0.22
set cartoon_oval_length, 1.18
set cartoon_oval_width, 0.30
set cartoon_highlight_color, grey85
set transparency_mode, 1
set_color query_blue, [0.141176, 0.454902, 0.823529]
set_color target_red, [0.941176, 0.356863, 0.384314]

# d1okra_d1m98a1: store one shared view for Full and w/o V.
load /home/ziting/Code/StructureAlignment/EAlign_clean/test_results/02_score_ablation/mechanism_explore/extrem_cases/d1okra_d1m98a1/full_ealign/d1okra_d1m98a1.pdb, okra_full
load /home/ziting/Code/StructureAlignment/EAlign_clean/test_results/02_score_ablation/mechanism_explore/extrem_cases/d1okra_d1m98a1/without_vector/d1okra_d1m98a1.pdb, okra_wo_v
hide everything, all
show cartoon, okra_full or okra_wo_v
color query_blue, (okra_full or okra_wo_v) and chain A
color target_red, (okra_full or okra_wo_v) and chain B
set cartoon_transparency, 0.0, (okra_full or okra_wo_v) and chain A
set cartoon_transparency, 0.08, (okra_full or okra_wo_v) and chain B
orient (okra_full or okra_wo_v) and chain B
zoom okra_full or okra_wo_v, 5.0
view okra_shared, store

disable okra_wo_v
view okra_shared, recall
set suspend_updates, 0
ray 1800, 1400
png /home/ziting/Code/StructureAlignment/EAlign_clean/test_results/02_score_ablation/mechanism_explore/extrem_cases/d1okra_d1m98a1/superposition/d1okra_d1m98a1_full_ealign.png, dpi=300

set suspend_updates, 1
disable okra_full
enable okra_wo_v
view okra_shared, recall
set suspend_updates, 0
ray 1800, 1400
png /home/ziting/Code/StructureAlignment/EAlign_clean/test_results/02_score_ablation/mechanism_explore/extrem_cases/d1okra_d1m98a1/superposition/d1okra_d1m98a1_without_vector.png, dpi=300

# d1is8a_d1umwa1: store one shared view for Full and w/o 3Di.
reinitialize
set quiet, 1
set suspend_updates, 1
bg_color white
set ray_opaque_background, 0
set ray_shadows, 0
set orthoscopic, 1
set antialias, 2
set ray_trace_mode, 1
set ambient, 0.46
set direct, 0.58
set specular, 0.28
set shininess, 24
set depth_cue, 0
set cartoon_fancy_helices, 1
set cartoon_smooth_loops, 1
set cartoon_flat_sheets, 1
set cartoon_sampling, 14
set cartoon_loop_radius, 0.22
set cartoon_oval_length, 1.18
set cartoon_oval_width, 0.30
set cartoon_highlight_color, grey85
set transparency_mode, 1
set_color query_blue, [0.141176, 0.454902, 0.823529]
set_color target_red, [0.941176, 0.356863, 0.384314]

load /home/ziting/Code/StructureAlignment/EAlign_clean/test_results/02_score_ablation/mechanism_explore/extrem_cases/d1is8a_d1umwa1/full_ealign/d1is8a_d1umwa1.pdb, is8_full
load /home/ziting/Code/StructureAlignment/EAlign_clean/test_results/02_score_ablation/mechanism_explore/extrem_cases/d1is8a_d1umwa1/without_3di/d1is8a_d1umwa1.pdb, is8_wo_3di
hide everything, all
show cartoon, is8_full or is8_wo_3di
color query_blue, (is8_full or is8_wo_3di) and chain A
color target_red, (is8_full or is8_wo_3di) and chain B
set cartoon_transparency, 0.0, (is8_full or is8_wo_3di) and chain A
set cartoon_transparency, 0.08, (is8_full or is8_wo_3di) and chain B
orient (is8_full or is8_wo_3di) and chain B
zoom is8_full or is8_wo_3di, 5.0
view is8_shared, store

disable is8_wo_3di
view is8_shared, recall
set suspend_updates, 0
ray 1800, 1400
png /home/ziting/Code/StructureAlignment/EAlign_clean/test_results/02_score_ablation/mechanism_explore/extrem_cases/d1is8a_d1umwa1/superposition/d1is8a_d1umwa1_full_ealign.png, dpi=300

set suspend_updates, 1
disable is8_full
enable is8_wo_3di
view is8_shared, recall
set suspend_updates, 0
ray 1800, 1400
png /home/ziting/Code/StructureAlignment/EAlign_clean/test_results/02_score_ablation/mechanism_explore/extrem_cases/d1is8a_d1umwa1/superposition/d1is8a_d1umwa1_without_3di.png, dpi=300

quit
