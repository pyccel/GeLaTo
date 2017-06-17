{ 
	{ node = "identifier", name = "Omega",
		 typeof = "Domain",
		 attributs = { dim = 2, kind = 'structured'}
 		},
	{ node = "identifier", name = "V",
		 typeof = "Space",
		 attributs = { domain = "Omega", kind = 'h1'}
 		},
	{ node = "identifier", name = "u",
		 typeof = "Field",
		 attributs = { space = "V"}
 		},
	{ node = "identifier", name = "a",
		 typeof = "BilinearForm",
		 attributs = { trial_space = "V", test_space = "V", fields = {}, kernel_name = "build_matrix_a" , filename = "kernels.lua",
		 n_rows = 1, n_cols = 1},
		 glt_symbol = "{ { {1},'*',{s1},'*',{m2} },'+',{ {1},'*',{m1},'*',{s2} } }"
 		},
	{ node = "identifier", name = "b",
		 typeof = "LinearForm",
		 attributs = { space = "V", fields = {}, kernel_name = "build_vector_b" , filename = "kernels.lua",
		 n_rows = 1, n_cols = 1}
 		}
 }