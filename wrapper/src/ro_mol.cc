#include "rust/cxx.h"
#include <GraphMol/GraphMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <DataStructs/ExplicitBitVect.h>
#include <GraphMol/Fingerprints/Fingerprints.h>
#include <GraphMol/MolStandardize/Tautomer.h>

namespace RDKit {
    using ExplicitBitVect = ::ExplicitBitVect;

    std::shared_ptr<ROMol> copy_mol(std::shared_ptr<ROMol> mol) {
        return std::shared_ptr<ROMol>(new ROMol(*mol));
    }

    std::shared_ptr<ROMol> smiles_to_mol(const std::string &smiles) {
        ROMol *mol = SmilesToMol(smiles);

        return std::shared_ptr<ROMol>(mol);
    }

     rust::String mol_to_smiles(std::shared_ptr<ROMol> mol) {
        return MolToSmiles(*mol);
    }
}