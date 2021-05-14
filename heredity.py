import csv
import itertools
import sys

PROBS = {

    # Unconditional probabilities for having gene
    "gene": {
        2: 0.01,
        1: 0.03,
        0: 0.96
    },

    "trait": {

        # Probability of trait given two copies of gene
        2: {
            True: 0.65,
            False: 0.35
        },

        # Probability of trait given one copy of gene
        1: {
            True: 0.56,
            False: 0.44
        },

        # Probability of trait given no gene
        0: {
            True: 0.01,
            False: 0.99
        }
    },

    # Mutation probability
    "mutation": 0.01
}


def main():

    # Check for proper usage
    if len(sys.argv) != 2:
        sys.exit("Usage: python heredity.py data.csv")
    people = load_data(sys.argv[1])

    # Keep track of gene and trait probabilities for each person
    probabilities = {
        person: {
            "gene": {
                2: 0,
                1: 0,
                0: 0
            },
            "trait": {
                True: 0,
                False: 0
            }
        }
        for person in people
    }

    # Loop over all sets of people who might have the trait
    names = set(people)
    for have_trait in powerset(names):

        # Check if current set of people violates known information
        fails_evidence = any(
            (people[person]["trait"] is not None and
             people[person]["trait"] != (person in have_trait))
            for person in names
        )
        if fails_evidence:
            continue

        # Loop over all sets of people who might have the gene
        for one_gene in powerset(names):
            for two_genes in powerset(names - one_gene):

                # Update probabilities with new joint probability
                p = joint_probability(people, one_gene, two_genes, have_trait)
                update(probabilities, one_gene, two_genes, have_trait, p)

    # Ensure probabilities sum to 1
    normalize(probabilities)

    # Print results
    for person in people:
        print(f"{person}:")
        for field in probabilities[person]:
            print(f"  {field.capitalize()}:")
            for value in probabilities[person][field]:
                p = probabilities[person][field][value]
                print(f"    {value}: {p:.4f}")


def load_data(filename):
    """
    Load gene and trait data from a file into a dictionary.
    File assumed to be a CSV containing fields name, mother, father, trait.
    mother, father must both be blank, or both be valid names in the CSV.
    trait should be 0 or 1 if trait is known, blank otherwise.
    """
    data = dict()
    with open(filename) as f:
        reader = csv.DictReader(f)
        for row in reader:
            name = row["name"]
            data[name] = {
                "name": name,
                "mother": row["mother"] or None,
                "father": row["father"] or None,
                "trait": (True if row["trait"] == "1" else
                          False if row["trait"] == "0" else None)
            }
    return data


def powerset(s):
    """
    Return a list of all possible subsets of set s.
    """
    s = list(s)
    return [
        set(s) for s in itertools.chain.from_iterable(
            itertools.combinations(s, r) for r in range(len(s) + 1)
        )
    ]


def joint_probability(people, one_gene, two_genes, have_trait):
    """
    Compute and return a joint probability.

    The probability returned should be the probability that
        * everyone in set `one_gene` has one copy of the gene, and
        * everyone in set `two_genes` has two copies of the gene, and
        * everyone not in `one_gene` or `two_gene` does not have the gene, and
        * everyone in set `have_trait` has the trait, and
        * everyone not in set` have_trait` does not have the trait.
    """
    #initialising probabilities that each people has a trait
    people_proba = dict.fromkeys(people.keys(), 0)
    people_trait = {}
    people_genes = dict.fromkeys(people.keys(), 0)
    total_proba_data = {}
    for each in people:
        total_proba_data[each] = {
            0: {
                True:None,
                False:None
            },
            1: {
                True: None,
                False:None
            },
            2: {
                True: None,
                False:None
            }
        }
    #initialising trait dictionary by setting people to True or False depending on they have the trait or not
    for each in people:
        have = False
        if each in have_trait:
            have = True
        people_trait[each] = have
    #update the number of genes possessed by each individual
    for each in one_gene:
        people_genes[each] = 1

    for each in two_genes:
        people_genes[each] = 2

    #iterate over all people and set the probability dictionary for people without parents
    for individual in people:
        #get the data of the individual
        individual_data = people[individual]

        # if the individual has no father and no mother
        if individual_data["father"] == None and individual_data["mother"] == None:

            trait = people_trait[individual]
            gene_count = people_genes[individual]
            proba = PROBS["trait"][gene_count][trait]

            unconditional_proba = PROBS["gene"][gene_count]
            people_proba[individual] = proba * unconditional_proba
            total_proba_data[individual][gene_count][trait] = proba * unconditional_proba
    #Now set the probability data for the people having parents
    for individual in people:
        individual_data = people[individual]
        # get the names of father and mother


        father = individual_data["father"]
        mother = individual_data["mother"]

        # if the individual has mother and father
        if father != None and mother != None:
            genes_count = people_genes[individual]
            trait = people_trait[individual]
            resultant_proba = 0
            # We can look at the trait from mother for 'gene_count_mom' copy of the gene
            gene_count_mom = people_genes[mother]
            mom_prob = 0

            # Let's calculate the probability given that he has '1' gene
            if gene_count_mom == 0:
                mom_prob = PROBS["mutation"]
            elif gene_count_mom == 2:
                mom_prob = 1-PROBS["mutation"]
            else:
                mom_prob = 0.5

            # Now we can look at the trait from mother 

            gene_count_dad = people_genes[father]
            dad_prob = 0

            if gene_count_dad == 0:
                dad_prob = PROBS["mutation"]
            elif gene_count_dad == 2:
                dad_prob = 1 - PROBS["mutation"]
            else:
                dad_prob = 0.5
            proba_gene = 0
            if genes_count==1:
                proba_gene = mom_prob * (1-dad_prob) + dad_prob * (1-mom_prob)
            elif genes_count == 2:
                proba_gene = mom_prob * dad_prob
            elif genes_count == 0:
                proba_gene = (1-mom_prob)*(1-dad_prob)
            prob_not_have_trait = PROBS["trait"][genes_count][trait]


            resultant_proba = proba_gene * prob_not_have_trait

            total_proba_data[individual][genes_count][trait] = resultant_proba
            people_proba[individual] = resultant_proba
            #get the probability for the probability of having the trait given he/she has 'genes_count' genes
    joint_proba = 1
    for h in people_proba:
        joint_proba *= people_proba[h]
    return joint_proba


def update(probabilities, one_gene, two_genes, have_trait, p):
    """
    Add to `probabilities` a new joint probability `p`.
    Each person should have their "gene" and "trait" distributions updated.
    Which value for each distribution is updated depends on whether
    the person is in `have_gene` and `have_trait`, respectively.
    """
    for individual in probabilities:
        individual_data = probabilities[individual]
        gene_number = 0
        if individual in one_gene:
            gene_number = 1
        elif individual in two_genes:
            gene_number = 2
        probabilities[individual]["gene"][gene_number] += p 
        if individual in have_trait:
            probabilities[individual]["trait"][True] += p
        elif individual not in have_trait:
            probabilities[individual]["trait"][False] += p
    #raise NotImplementedError


def normalize(probabilities):
    """
    Update `probabilities` such that each probability distribution
    is normalized (i.e., sums to 1, with relative proportions the same).
    """
    for individual in probabilities:
        trait_data = probabilities[individual]["trait"]
        gene_data = probabilities[individual]["gene"]
        # Normalizing trait data
        a = probabilities[individual]["trait"][True]
        b = probabilities[individual]["trait"][False]
        l = 1.0/(a+b)
        a = l*a
        b = l*b

        probabilities[individual]["trait"][True] = a
        probabilities[individual]["trait"][False]= b

        #Normalizing gene data
        a = gene_data[0]
        b = gene_data[1]
        c = gene_data[2]

        # normal data1
        k = 1.0/(a+b+c)
        a = k*a
        b = k*b
        c = k*c
        probabilities[individual]["gene"][0] = a
        probabilities[individual]["gene"][1] = b
        probabilities[individual]["gene"][2] = c


if __name__ == "__main__":
    main()