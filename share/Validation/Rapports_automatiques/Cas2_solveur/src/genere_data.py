dico_solveur={
    "gmres_3pts": " gmres { diag seuil -1 nb_it_max 5 }",
    "gmres_15pts": " gmres { diag seuil -1 nb_it_max 15 }",
    "gmres_15pts_KR10": " gmres { diag seuil -1 nb_it_max 15  dim_espace_krilov 10 }",
    "seuil_rel":" petsc bicgstab { seuil 1e-12 seuil_relatif 1e-3 impr precond  pilut { level 30 epsilon 0.01   } }",
    "seuil_rel_level20":" petsc bicgstab { seuil 1e-12 seuil_relatif 1e-3 impr precond  pilut { level 30 epsilon 0.01   } }",
    
    "seuil_non_rel":" petsc bicgstab { seuil 1e-12 seuil_relatif 0e-3 impr precond  pilut { level 30 epsilon 0.01   } }",
 #   "bicgstab_lapack":"gen { seuil 1.0e-15 solv_elem bicgstab precond ilu { type    2 filling 0 } impr  } ",
 #   "bicgstab_lapack_filling10":"gen { seuil 1.0e-15 solv_elem bicgstab precond ilu { type    2 filling 10 } impr } ",
 #   "bicgstab_lapack_gmres":"gen { seuil 1.0e-15 solv_elem gmres precond ilu { type    2 filling 10 } impr } ",

    }

if __name__=="__main__":
    import sys
    facsec=sys.argv[1]
    #solveur=sys.argv[3]

    f=open("3especes.data","r")
    lines=f.readlines()
    f.close()
    for solveur in dico_solveur.keys():
        
        f2=open("genere_"+solveur+"_"+facsec+".data","w")
        for line in lines:
            f2.write(line.replace("__solveur__",dico_solveur[solveur]).replace("__facsec__",facsec))
            pass
        f2.close()
        pass
