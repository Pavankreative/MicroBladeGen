from microbladegen.aero import Airfoil2D, plot2D_compare

# match with LPT vane 3, 50%
tauFact = 0.7
lengthFact = 42.0
axialChord = 0.040*lengthFact
radiusTE = 0.0004*lengthFact

# --- Top-to-bottom convention (default) ---
betaIn_ttb = -45.
betaOut_ttb = 72.
staggerAngle1 = 24.

vane1 = Airfoil2D(alpha1=betaIn_ttb, alpha2=betaOut_ttb, axial_chord=axialChord, stagger=staggerAngle1,left_to_right=True)
vane1.add_le_thickness(0.02)
vane1.add_ps_thickness(thicknessArray=[0.050*tauFact, 0.020*tauFact, -0.010*tauFact], expansion_ratio=1.2)
vane1.add_ss_thickness(thicknessArray=[0.150*tauFact, 0.400*tauFact, 0.340*tauFact, 0.160*tauFact], camberPercent=0.75, expansion_ratio=1.2)
vane1.match_le_thickness()
vane1.te_create(radius=radiusTE, wedge_ss=3.5, wedge_ps=2.4)
vane1.add_ss_flow_guidance_2(s_c=0.75, n=10)
vane1.center_le()

# --- Left-to-right convention ---
# LR angles measured from +X axis (axial), positive = above horizontal
# Equivalent to TTB(45, 72, 24): alpha1_lr = -45, alpha2_lr = 72
betaIn_lr = -45.
betaOut_lr = 72.
staggerAngle2 = 50.

vane2 = Airfoil2D(alpha1=betaIn_lr, alpha2=betaOut_lr, axial_chord=axialChord, stagger=staggerAngle2, left_to_right=True)
vane2.add_le_thickness(0.02)
vane2.add_ps_thickness(thicknessArray=[0.050*tauFact, 0.020*tauFact, -0.010*tauFact], expansion_ratio=1.2)
vane2.add_ss_thickness(thicknessArray=[0.150*tauFact, 0.400*tauFact, 0.340*tauFact, 0.160*tauFact], camberPercent=0.75, expansion_ratio=1.2)
vane2.match_le_thickness()
vane2.te_create(radius=radiusTE, wedge_ss=3.5, wedge_ps=2.4)
vane2.add_ss_flow_guidance_2(s_c=0.75, n=10)
vane2.center_le()

# Compare: TTB rotated 90° should match LR
vane2.center_le()
plot2D_compare([vane1, vane2], labels=['Vane 1', 'Vane 2'])
