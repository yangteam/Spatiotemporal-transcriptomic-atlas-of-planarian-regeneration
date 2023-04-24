spatial_colocalization_degree <- function(ligand, receptor, count_data){
    RL_int <- data.frame(ligand = count_data[ligand,],
                         receptor = count_data[receptor,])
    glmnb_rl <- glm.nb(receptor~ligand, data = RL_int)
    scd <- as.numeric(glmnb_rl$coefficients[2])
    return(scd)
}

