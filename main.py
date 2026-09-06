VALID_AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWY")

AMINO_ACID_ORDER = [
    "A", "C", "D", "E", "F", "G", "H", "I", "K", "L",
    "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"
]

amino_acid_weights = {
    "A": 89.09,
    "C": 121.15,
    "D": 133.10,
    "E": 147.13,
    "F": 165.19,
    "G": 75.07,
    "H": 155.16,
    "I": 131.17,
    "K": 146.19,
    "L": 131.17,
    "M": 149.21,
    "N": 132.12,
    "P": 115.13,
    "Q": 146.15,
    "R": 174.20,
    "S": 105.09,
    "T": 119.12,
    "V": 117.15,
    "W": 204.23,
    "Y": 181.19
}

amino_acid_pka_values = {
    "D": 3.9,
    "E": 4.1,
    "C": 8.3,
    "Y": 10.1,
    "H": 6.0,
    "K": 10.5,
    "R": 12.5
}

N_TERMINAL_PKA = 9.0
C_TERMINAL_PKA = 2.0

hydropathy_scale = {
    "A": 1.8,
    "C": 2.5,
    "D": -3.5,
    "E": -3.5,
    "F": 2.8,
    "G": -0.4,
    "H": -3.2,
    "I": 4.5,
    "K": -3.9,
    "L": 3.8,
    "M": 1.9,
    "N": -3.5,
    "P": -1.6,
    "Q": -3.5,
    "R": -4.5,
    "S": -0.8,
    "T": -0.7,
    "V": 4.2,
    "W": -0.9,
    "Y": -1.3
}

def validate_protein_sequence(sequence):
    sequence = sequence.replace(" ", "").upper()

    if not sequence:
        return False

    return all(amino_acid in VALID_AMINO_ACIDS for amino_acid in sequence)

def analyze_protein_sequence(sequence):
    sequence = sequence.replace(" ", "").upper()
    amino_acid_counts = {amino_acid: 0 for amino_acid in AMINO_ACID_ORDER}
    length = len(sequence)

    for amino_acid in sequence:
        if amino_acid in VALID_AMINO_ACIDS:
            amino_acid_counts[amino_acid] += 1

    return amino_acid_counts, length

def calculate_molecular_weight(sequence):
    sequence = sequence.replace(" ", "").upper()
    molecular_weight = 0.0

    for amino_acid in sequence:
        if amino_acid in amino_acid_weights:
            molecular_weight += amino_acid_weights[amino_acid]

    return molecular_weight

def estimate_isoelectric_point(sequence):
    sequence = sequence.replace(" ", "").upper()
    low_pH = 0.0
    high_pH = 14.0
    max_iterations = 1000

    for _ in range(max_iterations):
        pH = (low_pH + high_pH) / 2
        net_charge = 0.0

        # N-terminal contribution
        net_charge += 1 / (1 + 10 ** (pH - N_TERMINAL_PKA))

        # C-terminal contribution
        net_charge -= 1 / (1 + 10 ** (C_TERMINAL_PKA - pH))

        for amino_acid in sequence:
            if amino_acid in amino_acid_pka_values:
                pKa = amino_acid_pka_values[amino_acid]
                if amino_acid in "DECY":
                    net_charge -= 1 / (1 + 10 ** (pKa - pH))
                elif amino_acid in "HKR":
                    net_charge += 1 / (1 + 10 ** (pH - pKa))

        if abs(net_charge) < 0.01:
            break

        if net_charge > 0:
            low_pH = pH
        else:
            high_pH = pH

    return round(pH, 2)

def calculate_hydropathy_index(sequence):
    sequence = sequence.replace(" ", "").upper()
    total_hydropathy = 0.0
    length = len(sequence)

    for amino_acid in sequence:
        if amino_acid in hydropathy_scale:
            total_hydropathy += hydropathy_scale[amino_acid]

    if length == 0:
        return 0.0

    return total_hydropathy / length

def interpret_hydropathy_index(hydropathy_index):
    if hydropathy_index > 0:
        return "Hydrophobic"
    elif hydropathy_index < 0:
        return "Hydrophilic"
    else:
        return "Neutral"

def calculate_amino_acid_groups(sequence):
    sequence = sequence.replace(" ", "").upper()
    groups = {
        "Nonpolar": 0,
        "Polar": 0,
        "Acidic": 0,
        "Basic": 0
    }

    for amino_acid in sequence:
        if amino_acid in "ACFGILMPVWY":
            groups["Nonpolar"] += 1
        elif amino_acid in "DE":
            groups["Acidic"] += 1
        elif amino_acid in "HKR":
            groups["Basic"] += 1
        elif amino_acid in "NQST":
            groups["Polar"] += 1

    return groups

def interpret_charge_tendency(sequence):
    sequence = sequence.replace(" ", "").upper()
    acidic_count = sum(sequence.count(aa) for aa in "DE")
    basic_count = sum(sequence.count(aa) for aa in "HKR")

    if acidic_count > basic_count:
        return "Acidic"
    elif basic_count > acidic_count:
        return "Basic"
    else:
        return "Neutral"

def clean_protein_sequence(sequence):
    sequence = sequence.upper()

    lines = sequence.splitlines()

    if lines and lines[0].startswith(">"):
        lines = lines[1:]

    sequence = "".join(lines)
    sequence = "".join(sequence.split())

    return sequence

def main():
    print("Protein Sequence Analyzer")
    print("=========================")

    print("Enter a protein sequence: ")
    print("Please enter the sequence in FASTA format or as a plain string. The program will clean and validate the input.")

    sequence_lines = []

    while True:
        line = input()
        if line.strip() == "":
            break

        sequence_lines.append(line)

    sequence = "\n".join(sequence_lines)

    sequence = clean_protein_sequence(sequence)

    if validate_protein_sequence(sequence):
        print("Valid protein sequence!")

        amino_acid_counts, length = analyze_protein_sequence(sequence)
        molecular_weight = calculate_molecular_weight(sequence)
        isoelectric_point = estimate_isoelectric_point(sequence)
        hydropathy_index = calculate_hydropathy_index(sequence)
        hydropathy_interpretation = interpret_hydropathy_index(hydropathy_index)
        amino_acid_groups = calculate_amino_acid_groups(sequence)
        charge_tendency = interpret_charge_tendency(sequence)
        print(f"Length of the sequence: {length}")
        print(f"Molecular Weight: {molecular_weight:.2f} Da")
        print(f"Isoelectric Point: {isoelectric_point}")
        print(f"Hydropathy Index: {hydropathy_index:.2f}")
        print(f"Hydropathy Interpretation: {hydropathy_interpretation}")
        print(f"Charge Tendency: {charge_tendency}")

        print("\nAmino Acid Counts:")
        for amino_acid, count in amino_acid_counts.items():
            if count > 0:
                percentage = (count / length) * 100
                print(f"{amino_acid}: {count} ({percentage:.2f}%)")

        print("\nAmino Acid Groups:")        
        for group, count in amino_acid_groups.items():
            if count > 0:
                percentage = (count / length) * 100
                print(f"{group}: {count} ({percentage:.2f}%)")
    else:
        print("Invalid protein sequence.")


if __name__ == "__main__":
    main()

