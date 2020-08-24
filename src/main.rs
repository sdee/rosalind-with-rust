use counter::Counter;
use regex::Regex;
use std::collections::HashMap;
use std::collections::HashSet;
#[macro_use]
extern crate maplit;

// Problem BA1E: find Clumps in String
pub fn find_clumps() {
    let s = "CGGACTCGACAGATGTGAAGAAATGTGAAGACTGAGTGAAGAGAAGAGGAAACACGACACGACATTGCGACATAATGTACGAATGTAATGTGCCTATGGC";
    let chars: Vec<char> = s.chars().collect();
    let counter = chars
        .windows(5)
        .map(|w| w.into_iter().collect::<String>())
        .collect::<Counter<_>>();
    let clumps = counter
        .iter()
        .filter(|&(key, value)| *value > 3)
        .collect::<Vec<_>>();
    for (key, val) in clumps.iter() {
        println!("{:?}", key);
    }
}

pub fn generate_d_neighborhood() {
    //
}

pub fn transcribe(s: &str) -> std::string::String {
    return s.replace("T", "U");
}

// Problem: rvec
pub fn reverse_complement(s: &str) -> std::string::String {
    let complements = hashmap! {'A'=>'T', 'T'=>'A', 'C'=>'G', 'G'=>'C'};
    let reverse_complement = s.chars().rev().map(|c| complements[&c]).collect::<String>();
    return reverse_complement;
}

pub fn translate(s: &str) {
    let codons = hashmap! {"UUU"=>"F", "UUC"=>"F", "UUA"=>"L", "UUG"=>"L", "UCU"=>"S", "UCC"=>"S", "UCA"=>"S", "UCG"=>"S", "UAU"=>"Y", "UAC"=>"Y", "UGU"=>"C", "UGC"=>"C", "UGG"=>"W", "CUU"=>"L", "CUC"=>"L", "CUA"=>"L", "CUG"=>"L", "CCU"=>"P", "CCC"=>"P", "CCA"=>"P", "CCG"=>"P", "CAU"=>"H", "CAC"=>"H", "CAA"=>"Q", "CAG"=>"Q", "CGU"=>"R", "CGC"=>"R", "CGA"=>"R", "CGG"=>"R", "AUU"=>"I", "AUC"=>"I", "AUA"=>"I", "AUG"=>"M", "ACU"=>"T", "ACC"=>"T", "ACA"=>"T", "ACG"=>"T", "AAU"=>"N", "AAC"=>"N", "AAA"=>"K", "AAG"=>"K", "AGU"=>"S", "AGC"=>"S", "AGA"=>"R", "AGG"=>"R", "GUU"=>"V", "GUC"=>"V", "GUA"=>"V", "GUG"=>"V", "GCU"=>"A", "GCC"=>"A", "GCA"=>"A", "GCG"=>"A", "GAU"=>"D", "GAC"=>"D", "GAA"=>"E", "GAG"=>"E", "GGU"=>"G", "GGC"=>"G", "GGA"=>"G", "GGG"=>"G"};

    let start_codon = "AUG";
    let stop_codons = ["UAA", "UAG", "UGA"];
    let mut protein_strings = Vec::new();
    for window in s.chars().collect::<Vec<char>>().chunks(3) {
        let codon = window
            .into_iter()
            .map(|i| i.to_string())
            .collect::<Vec<_>>()
            .join("");
        if stop_codons.contains(&&codon.as_str()) {
            break;
        } else if codon.as_str() == start_codon {
            &protein_strings.push("".to_string());
        }
        let aa = &codons.get(codon.as_str()).unwrap_or(&"").to_string();
        for protein_string in protein_strings.iter_mut() {
            *protein_string = protein_string.to_owned() + &aa;
        }
    }
    if protein_strings.len() > 0 {
        for protein_string in protein_strings {
            println!("{:?}", protein_string);
        }
    }
}

fn orfs(s: &str) {
    let rna = transcribe(&s);
    let rev = reverse_complement(&s);
    let rev_rna = transcribe(&rev);
    let ORFs = [
        &rev_rna,
        &rev_rna[1..],
        &rev_rna[2..],
        &rna,
        &rna[1..],
        &rna[2..],
    ];
    for ORF in ORFs.iter() {
        translate(ORF);
    }
}

// Problem ID: DNA
pub fn count_nucleotides(s: &String) {
    let counter = s.chars().collect::<Counter<_>>();
    let base_count = counter.into_map();
    let a_count = if base_count.contains_key(&'A') {
        base_count[&'A']
    } else {
        0
    };
    let c_count = if base_count.contains_key(&'C') {
        base_count[&'C']
    } else {
        0
    };
    let g_count = if base_count.contains_key(&'G') {
        base_count[&'G']
    } else {
        0
    };
    let t_count = if base_count.contains_key(&'T') {
        base_count[&'T']
    } else {
        0
    };
    println!("{} {} {} {}", a_count, c_count, g_count, t_count);
}

// Problem ID: HAMM
pub fn calculate_hamming_distance(s1: &String, s2: &String) {
    let mut hamming_distance = 0;
    let iter = s1.chars().zip(s2.chars());
    for (i, j) in iter {
        if i != j {
            hamming_distance += 1;
        }
    }
    println!("{}", hamming_distance);
}

// Problem ID: GC
pub fn calculate_gc_content(text: &String) {
    let lines = text.lines();
    let mut gc_by_sequence = HashMap::new();
    let mut identifier = "";
    for line in lines {
        let rosalind_re = Regex::new(r".*>(.*)").unwrap();
        if line.contains(">") {
            match rosalind_re.captures(line) {
                Some(x) => identifier = x.at(1).unwrap(),
                None => unreachable!(),
            }
        } else {
            let base_count = line.chars().collect::<Counter<_>>();
            let gc_content = base_count[&'G'] + base_count[&'C'];
            let total_num_bases = gc_content + base_count[&'A'] + base_count[&'T'];
            let gc_score: f64 = gc_content as f64 / total_num_bases as f64;
            gc_by_sequence.insert(identifier, gc_score);
        }
    }
    let highest_gc_sequence = gc_by_sequence
        .iter()
        .max_by(|a, b| a.1.partial_cmp(b.1).unwrap());
    println!(
        "{:?} has highest GC content with a GC content {}",
        highest_gc_sequence.unwrap().0,
        highest_gc_sequence.unwrap().1
    );
}

// Problem ID: TRAN
pub fn calculate_transition_and_transversion(s1: &String, s2: &String) {
    let purines: HashSet<char> = "AG".chars().map(|c| c).collect();
    let pyrimidines: HashSet<char> = "CT".chars().map(|c| c).collect();
    let mut number_transitions = 0;
    let mut number_transversions = 0;
    let iter = s1.chars().zip(s2.chars());
    for (i, j) in iter {
        if i != j {
            if (purines.contains(&i) && purines.contains(&j))
                || (pyrimidines.contains(&i) && pyrimidines.contains(&j))
            {
                number_transitions += 1;
            } else {
                number_transversions += 1;
            }
        }
    }
    let transition_transversion_ratio: f64 =
        number_transitions as f64 / number_transversions as f64;
    println!(
        "Transition/transversion ratio: {:?}",
        transition_transversion_ratio
    );
}

// Problem: SSEQ
pub fn find_splice_motif(s1: &String, s2: &String) {
    let mut motif_iter = s2.chars();
    let mut char_need_to_match = motif_iter.next();
    let mut splice_coords: Vec<String> = Vec::new();
    for (i, c) in s1.chars().enumerate() {
        if c == char_need_to_match.unwrap() {
            splice_coords.push(i.to_string());
            if motif_iter.as_str() == "" {
                println!("Coords: {}", splice_coords.join(" "));
                break;
            }
            char_need_to_match = motif_iter.next();
        }
    }
}

// Problem: Generate the k-mer Composition of a String
pub fn generate_kmer_composition(s1: &String, k: usize) {
    println!("KMER");
    println!("{}", s1);
    let chars: Vec<char> = s1.chars().collect();
    for window in chars.windows(k) {
        let s: String = window.into_iter().collect();
        println!("{}", s);
        // let s2: String = window.into_iter().map(|i| i.to_string()).collect::<Vec<_>>().join(" ");
        // println!("{}", s2);
    }
}

// Problem: Reconstruct a String from its Genome Path
pub fn reconstruct_sequence_from_kmer() {
    let kmers = [
        "ACCGA".to_string(),
        "CCGAA".to_string(),
        "CGAAG".to_string(),
        "GAAGC".to_string(),
        "AAGCT".to_string(),
    ];
    let mut sequence = String::new();
    let k: usize = kmers[0].len();
    for (i, kmer) in kmers.iter().enumerate() {
        if i == 0 {
            sequence += kmer;
        } else {
            sequence += &kmer.chars().nth(k - 1).unwrap().to_string();
        }
    }
}

// Given a collection of kmers, construct an adjacency list representing the overlap graph. Handles graphs that aren't a strict linked list but not cycles
pub fn construct_overlap_graph() {
    let kmers = [
        "ATGCG".to_string(),
        "GCATG".to_string(),
        "CATGC".to_string(),
        "AGGCA".to_string(),
        "GGCAT".to_string(),
    ];
    let mut suffixes = HashMap::new();
    let mut prefixes = HashMap::new();

    let sequence = String::new();
    let k: usize = kmers[0].len();

    for kmer in kmers.iter() {
        prefixes
            .entry(&kmer[0..k - 1])
            .or_insert(Vec::new())
            .push(kmer);
        suffixes.entry(&kmer[1..k]).or_insert(Vec::new()).push(kmer);
    }
    for (suffix, sequences_with_suffix) in &suffixes {
        let matches = &prefixes.get(suffix);
        for seq1 in sequences_with_suffix {
            if matches.is_some() {
                for seq2 in matches.unwrap().iter() {
                    println!("{} -> {}", seq1, seq2);
                }
            }
        }
    }
}

// generate kmers, create adjacency graph, print
pub fn construct_de_bruijn_graph_from_string() {
    let k: usize = 4;
    let sequence = "AAGATTCTCTAC".to_string();
    let mut prev_window: String = sequence[0..k - 1].to_string();
    let mut current_window: String;
    let mut graph = HashMap::new();
    let chars: Vec<char> = sequence.chars().collect();
    for window in chars[1..].windows(k - 1) {
        current_window = window.into_iter().collect();
        graph
            .entry(prev_window)
            .or_insert(Vec::new())
            .push(current_window.clone());
        prev_window = current_window;
    }
    for (seq1, value) in &graph {
        println!(
            "{}->{}",
            seq1,
            value
                .into_iter()
                .map(|i| i.to_string())
                .collect::<Vec<_>>()
                .join(", ")
        );
    }
}

fn main() {
    let s = String::from("AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC");
    count_nucleotides(&s);

    let s1 = String::from("GAGCCTACTAACGGGAT");
    let s2 = String::from("CATCGTAATGACGGCCT");
    calculate_hamming_distance(&s1, &s2);

    let text = String::from(
        ">Rosalind_6404
    CCTGCGGAAGATCGGCACTAGAATAGCCAGAACCGTTTCTCTGAGGCTTCCGGCCTTCCCTCCCACTAATAATTCTGAGG
    >Rosalind_5959
    CCATCGGTAGCGCATCCTTAGTCCAATTAAGTCCCTATCCAGGCGCTCCGCCGAAGGTCTATATCCATTTGTCAGCAGACACGC
    >Rosalind_0808
    CCACCCTCGTGGTATGGCTAGGCATTCAGGAACCGGAGAACGCTTCAGACCAGCCCGGACTGGGAACCTGCGGGCAGTAGGTGGAAT",
    );

    calculate_gc_content(&text);

    let s3 = String::from(
        "GCAACGCACAACGAAAACCCTTAGGGACTGGATTATTTCGTGATCGTTGTAGTTATTGGAAGTACGGGCATCAACCCAGTT",
    );
    let s4 = String::from(
        "TTATCTGACAAAGAAAGCCGTCAACGGCTGGATAATTTCGCGATCGTGCTGGTTACTGGCGGTACGAGTGTTCCTTTGGGT",
    );
    calculate_transition_and_transversion(&s3, &s4);

    let s5 = String::from("ACGTACGTGACG");
    let s6 = String::from("GTA");
    find_splice_motif(&s5, &s6);

    let s7 = String::from("ACGTCCCCACGTGACG");
    let s8 = String::from("GTA");
    find_splice_motif(&s7, &s8);

    let s9 = String::from("CAATCCAAC");
    generate_kmer_composition(&s9, 5);

    // let kmers = [
    //     "ACCGA".to_string(),
    //     "CCGAA".to_string(),
    //     "CGAAG".to_string(),
    //     "GAAGC".to_string(),
    //     "AAGCT".to_string(),
    // ];

    reconstruct_sequence_from_kmer();

    construct_overlap_graph();

    construct_de_bruijn_graph_from_string();

    let s10 = "AAAACCCGGT";
    reverse_complement(s10);
    println!("{:?}", transcribe("GATGGAACTTGACTACGTAAATT"));
    translate("AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA");
    orfs("AGCCATGTAGCTAACTCAGGTTACATGGGGATGACCCCGCGACTTGGATTAGAGTCTCTTTTGGAATAAGCCTGAATGATCCGAGTAGCATCTCAG");

    println!("{}", "-------");
    generate_kmer_composition(&String::from("ACGT"), 2);
    find_clumps()
}
