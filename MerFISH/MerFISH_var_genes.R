# Curated genes used in MerFISH


var_genes = c(
"Tmem119",
"Fcrls",
"Eif3f",
"Itm2b",
"Arhgdib",
"Selenop",
"Fth1",
"Cd81",
"Cd34",
"P2ry13",
"Zfhx3",
"Mertk",
"Irf2bp2",
"Pag1",
"Tnrc6b",
"Slco2b1",
"Pld1",
"Cx3cr1",
"Cd164",
"Oasl2",
"Rtp4",
"Abcg1",
"Spp1",
"Cd63",
"Ctsb",
"C3ar1",
"Apoe",
"Lyz2",
"Ccr1",
"Tmem176a",
"Tmem176b",
"Mcm4",
"Ube2c",
"Mki67",
"Stard8",
"Pou3f2",
"Htr7",
"Cartpt",
"Qrfpr",
"Ndst4",
"Rorb",
"Smoc2",
"Kcnh5",
"Dock5",
"Susd5",
"Shroom3",
"Tcerg1l",
"Fezf2",
"Bcl11b",
"Foxo1",
"Deptor",
"Trhr",
"Dkk2",
"Lgr5",
"Crym",
"Tle4",
"Foxp2",
"Cpa6",
"Sulf1",
"Ccn2",
"Moxd1",
"Nxph4",
"Slc17a7",
"Neurod2",
"Aldh1l1",
"Slc6a20a",
"Slc6a13",
"Col4a6",
"Sox10",
"Itpr2",
"Opalin",
"Ptgds",
"Anln",
"Btg2",
"Egr2",
"Ier2",
"Dusp1",
"Il33",
"Npsr1",
"Abca8a",
"Erbb4",
"Gabbr2",
"F13a1",
"Lyve1",
"Ccr2"
)

Microglia_genes=c(
"Tmem119",  # general
"Fcrls",
"Eif3f",  # hom1
"Itm2b",
"Arhgdib",
"Selenop",
"Fth1",
"Cd81",
"Cd34",
"P2ry13",  # hom2
"Zfhx3",
"Mertk",
"Irf2bp2",
"Pag1",
"Tnrc6b",
"Slco2b1",
"Pld1",
"Cx3cr1",
"Cd164",
"Oasl2",  # innate immune
"Rtp4",
"Abcg1",  # inflammatory
"Spp1",
"Cd63",
"Ctsb",
"C3ar1",
"Apoe",  # apoe
"Lyz2",
"Ccr1",  # ccr1
"Tmem176a",
"Tmem176b"
)

ExN_genes=c(
"Stard8",  # L2-3 CPN
"Pou3f2",
"Htr7",
"Cartpt",
"Qrfpr",
"Ndst4",
"Rorb",  # L4 Stellate
"Smoc2",
"Kcnh5",
"Dock5",  # L5-6 CPN
"Susd5",
"Shroom3",
"Tcerg1l",  # SCPN
"Fezf2",
"Bcl11b",
"Foxo1",
"Deptor",
"Trhr",
"Dkk2",
"Lgr5",
"Crym",
"Tle4",  # CThPN
"Foxp2",
"Cpa6",
"Sulf1",
"Ccn2",  # L6b
"Moxd1",
"Nxph4",
"Slc17a7",  # all
"Neurod2"
)

Endo_genes=c(
"Slc6a20a",
"Slc6a13",
"Col4a6"
)

Oligo_genes=c(
"Sox10",
"Itpr2",
"Opalin",
"Ptgds",
"Anln",
"Btg2",
"Egr2",
"Ier2",
"Dusp1",
"Il33",
"Npsr1",
"Abca8a"
)

Neuronal_genes = c(
"Stard8",  # L2-3 CPN
"Pou3f2",
"Htr7",
"Cartpt",
"Qrfpr",
"Ndst4",
"Rorb",  # L4 Stellate
"Smoc2",
"Kcnh5",
"Dock5",  # L5-6 CPN
"Susd5",
"Shroom3",
"Tcerg1l",  # SCPN
"Fezf2",
"Bcl11b",
"Foxo1",
"Deptor",
"Trhr",
"Dkk2",
"Lgr5",
"Crym",
"Tle4",  # CThPN
"Foxp2",
"Cpa6",
"Sulf1",
"Ccn2",  # L6b
"Moxd1",
"Nxph4",
"Slc17a7",  # all
"Neurod2",
"Erbb4",  # IN
"Gabbr2"
)

NonNeuronal_genes=c(
"Tmem119",  # general
"Fcrls",
"Eif3f",  # hom1
"Itm2b",
"Arhgdib",
"Selenop",
"Fth1",
"Cd81",
"Cd34",
"P2ry13",  # hom2
"Zfhx3",
"Mertk",
"Irf2bp2",
"Pag1",
"Tnrc6b",
"Slco2b1",
"Pld1",
"Cx3cr1",
"Cd164",
"Oasl2",  # innate immune
"Rtp4",
"Abcg1",  # inflammatory
"Spp1",
"Cd63",
"Ctsb",
"C3ar1",
"Apoe",  # apoe
"Lyz2",
"Ccr1",  # ccr1
"Tmem176a",
"Tmem176b",
"Slc6a20a",  # Endothelial
"Slc6a13",
"Col4a6",
"Sox10",  # Oligo
"Itpr2",
"Opalin",
"Ptgds",
"Anln",
"Btg2",
"Egr2",
"Ier2",
"Dusp1",
"Il33",
"Npsr1",
"Abca8a",
"Aldh1l1"  # Astro
)

mg.target=c(
"P2ry6",
"Nrp1",
"Plxna4",
"Plxnb2",
"Sort1",
"Tgfbr2",
"Tgfbr1",
"Bmpr2",
"Bmpr1a",
"Bmpr1b",
"Nrp2",
"Epor",
"Nrp1",
"Epha2"
)

n.target=c(
"Wnt4",
"Nampt",
"Sema3a",
"Sema4g",
"Tnc",
"Bdnf",
"Tgfb2",
"Sema3c",
"Sema4c",
"Vegfa",
"Vegfb",
"Gas6",
"Epo",
"Efnb1"
)
