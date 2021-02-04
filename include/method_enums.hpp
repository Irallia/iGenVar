#pragma once

//!\brief An enum for the different clustering methods.
enum clustering_methods
{
    simple_clustering = 0,
    hierarchical_clustering = 1,
    self_balancing_binary_tree = 2,
    candidate_selection_based_on_voting = 3
};

//!\brief An enum for the different refinement methods.
enum refinement_methods
{
    no_refinement = 0,
    sViper_refinement_method = 1,
    sVirl_refinement_method = 2
};