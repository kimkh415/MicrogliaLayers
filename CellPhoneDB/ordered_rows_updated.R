# ordering and grouping of the interacting pairs found to be significant in at least one result

general = c(
"ACVR_1B2B receptor_GDF1",       # all
"TNFSF12_TNFRSF25",
"NOTCH1_JAG2",
"CLCF1_CNTF-1R",
"COL15A1_a11b1 complex",
"COL27A1_a11b1 complex",
"CCL8_ACKR1",
"EFNB1_EPHB3",
"IL16_GRIN2D",
"NOTCH2_JAG2",
"EFNB1_EPHB6",
"ADRB2_VEGFB",
"CKLF_LRP6",
"CCL23_IDE",
"JAM2_JAM3",
"CD83_FGFR2",
"GAS6_TYRO3",
"PROS1_TYRO3",
"ENTPD1_ADORA1",
"TGFB1_TGFBR3",
"a6b1 complex_LAMC1",
"PLXNB2_SEMA4D",
"P2RY6_COPA",
"PLXNB2_PTN",
"LAMP1_FAM3C",
"GRN_SORT1",
"IL16_GRIN2B",
"IL16_GRIN2A",
"CSF1R_IL34"
)

ko_specific = c(
"PLXNB2_SEMA4C",
"TGFB1_aVb6 complex",
"SORT1_BDNF",
"CRLF2_TSLPR"
)

L14 = c(
"NOTCH1_WNT4",
"EFNB1_EPHA4",
"ACVR_1B2B receptor_INHBA",
"ACVR_1A2B receptor_INHBA",
"PlexinA4_complex1_SEMA3A",
"PlexinA2_complex1_SEMA3A",
"NRP1_SEMA3A",
"PDGFB_ADGRV1",
"NECTIN2_NECTIN3",
"JAM2_a4b1 complex",
"a9b1 complex_TNC",
"XPR1_FGFR2",
"CADM1_CADM1"
)


mg.driven = c(
"TSLPR_CRLF2",
"NRP1_VEGFA",
"NRP1_VEGFB",
"EPHA2_EFNA5",
"MERTK_GAS6",
"LRP1_PDGFB",
"CX3CR1_CX3CL1"
)

L5 = c(
"TGFbeta receptor1_TGFB3",
"FZD4_WNT2B",
"FZD4_WNT2",
"IL1A_IL1 receptor inhibitor",
"TGFbeta receptor2_TGFB3",
"GAS6_MERTK",
"TGFB1_TGFbeta receptor1",
"BMR1B_AVR2B_BMP6",
"BMR1B_AVR2A_BMP6",
"BMR1A_AVR2B_BMP6",
"BMR1A_ACR2A_BMP6",
"BMPR1A_BMPR2_BMP6",
"ACVR_1A2B receptor_BMP6",
"ACVR_1A2A receptor_BMP6",
"TGFB1_TGFbeta receptor2",
"FZD4_WNT7B",
"IL15_IL15RA",
"ACVR_1B2B receptor_INHBB",
"ACVR_1A2B receptor_INHBB",
"PLXNB2_SEMA4G",
"TGFbeta receptor1_TGFB2",
"IL1A_IL1 receptor",
"PECAM1_CD38",
"LRP1_MDK",
"BMPR1A_BMPR1B_BMPR2",
"TGFbeta receptor2_TGFB2"
)

L6 = c(
"NRP2_SEMA3F",
"EFNB1_EPHB1",
"IL13 receptor_IL4",
"EPHA2_EFNA2",
"NECTIN2_CD226",
"EFNB1_EPHB2",
"FAM3C_GLRA2",
"VEGFB_NRP1",
"NRP2_SEMA3C"
)

L15 = c(
"SMO_WNT4",
"NRP2_VEGFA",
"NECTIN4_NECTIN1",
"F11R_BDNF",
"P2RY6_NAMPT",
"COPA_SORT1"
)
row.ord = list(General=general, KO=ko_specific, L14=L14, MgDriven=mg.driven, L5=L5, L6=L6, L15=L15)
