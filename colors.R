
# provides colors, but also preferred order for plotting cell types
ct_cols = function() {
    c(
        # Adult OL
        "OL (Adult)"='#FEB24C',
        # Adult GL
        "GL Myh1 (Adult)"='#E31A1C',
        "GL Myh4 (Adult)"='#800026',
        # Adult Slow
        "GL Slow (Adult)"='#FC4E2A',
        # Adult MTJ
        "MTJ (Adult)"='#3182BD',
        # Adult NMJ
        "NMJ (Adult)"='#D772D7',
        # Adult Olfml2a
        "Olfml2a (Adult)"='#636363',

        # P5
        "OL (P5)"='#8C76D7',
        "GL (P5)"='#BCBDDC',
        "Slow (P5)"='#807DBA',
        "MTJ (P5)"='#9ECAE1',
        "NMJ (P5)"='#C01ED7',
        "Olfml2a (P5)"='#969696',
        
        # E14
        "Pro-OL (E14)"='#00441B',
        "Pro-GL (E14)"='#238B45',
        "Pro-MTJ (E14)"='#41AB5D',
        "Pro-NMJ (E14)"='#C7E9C0',
        "Myoblasts (E14)"='#74C476'
        #"Myoblasts (E14)"='#74C476'
    )
}

stage_cols = function() {
    c(
        "Adult"='#E31A1C',
        "P5"='#BCBDDC',
        "E18"='#00441B',
        "E14"='#74C476'
    )
}

# provides colors, but also preferred order for plotting cell types
ct_cols_limb = function() {
    c(
        # Adult Myh1/Myh4
        "Myonuclei Myh1/Myh2 (Adult)"='#E31A1C',
        "Myonuclei Myh4 (Adult)"='#800026',
        # Adult Slow
        "Slow (Adult)"='#FC4E2A',
        # Adult NMJ
        "NMJ (Adult)"='#D772D7',
        # Adult MTJ
        "MTJ (Adult)"='#3182BD',
        
        # P5
        "Myonuclei Myh8 (P5)"='#BCBDDC',
        "Slow (P5)"='#807DBA',
        "NMJ (P5)"='#C01ED7',
        "MTJ (P5)"='#9ECAE1',
        
        
        # Myoblasts/E14/E18
        "Myonuclei Myh3/Myh8 (E18)"='#00441B',
        "Myonuclei Myh3/Myog (E18)"='#238B45',
        "Myonuclei Myh3 Col19a1 (E14)"='#C7E9C0',
        "Myonuclei Myh3 Col22a1 (E14)"='#41AB5D',
        
        "Myoblasts"='#74C476',
        "MuSC"='#8EAB0E'
    )
}

# provides colors, but also preferred order for plotting cell types
subtype_cols = function() {
    c(
        # Adult OL
        "Myh2 Pcdh15 (Adult)"='#FEB24C',
        "Myh15 (Adult)"='#FD8D3C',
        # Adult GL
        "Myh1/Myh2/Myh13 Mlycd (Adult)"='#E31A1C',
        "Myh1/Myh2/Myh13 Kcnq5 (Adult)"='#FED976',
        "Myh4 Mlycd (Adult)"='#800026',
        "Myh4 Kcnq5 (Adult)"='#BD0026',
        # Adult Slow
        "Myh7 (Adult)"='#FC4E2A',
        # Adult MTJ
        "MTJ Myh15 (Adult)"='#08519C',
        "MTJ (Adult)"='#3182BD',
        # Adult NMJ
        "NMJ Myh7/Myh8 (Adult)"='#D772D7',
        # Adult Olfml2a
        "Olfml2a (Adult)"='#636363',
        "Myh7/Olfml2a (Adult)"='#252525',

        # P5
        "Myh8 Myh2 (P5)"='#DADAEB',
        "Myh8 Myh2/Myh1 (P5)"='#BCBDDC',
        "Myh8 Myh4/Myh1 (P5)"='#9E9AC8',
        "Myh7 (P5)"='#807DBA',
        "Myh15 (P5)"='#6A51A3',
        "Myh8 (P5)"='#54278F',
        "Myh8 Myh1 (P5)"='#3F007D',
        "Myh8 Lmx1a/Myh2 (P5)"='#8C76D7',
        "MTJ (P5)"='#9ECAE1',
        "MTJ Myh15 (P5)"='#6BAED6',
        "NMJ Myh8 (P5)"='#C01ED7',
        "NMJ Myh7 (P5)"='#EB4CBF',
        "Olfml2a (P5)"='#969696',
        
        # E14
        "Myh3 Col25a1/Etv5/Cdh20 (NMJ-like, E14)"='#C7E9C0',
        "Myh3 Synpo2/Lpar1 (E14)"='#A1D99B',
        "Myh3 Sema6a/Adamts5 (E14)"='#74C476',
        "Myh3 Col22a1/Sox5/Chodl (MTJ-like, E14)"='#41AB5D',
        "Myh3 Col19a1/Tbx18/Gpc3 (E14)"='#238B45',
        "Myh3 Dcx/Col22a1/Scn3a (E14)"='#006D2C',
        "Myh3 Lmx1a/Hmga2/Stk32b (E14)"='#00441B',
        "Myh3 Myo16/Elmo1 (E14)"='#9FD645',
        "Myh3 Col19a1/Tll2/Vav3 (E14)"='#8EAB0E'
    )
}