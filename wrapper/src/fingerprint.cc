#include "rust/cxx.h"
#include <DataStructs/ExplicitBitVect.h>
#include <GraphMol/Fingerprints/Fingerprints.h>
#include <GraphMol/Fingerprints/MorganFingerprints.h>

namespace RDKit
{
    std::shared_ptr<ExplicitBitVect> fingerprint_mol(std::shared_ptr<ROMol> mol)
    {
        return std::shared_ptr<ExplicitBitVect>(RDKFingerprintMol(*mol));
    }

    std::unique_ptr<std::vector<uint64_t>> fingerprint_mol2(std::shared_ptr<ROMol> mol)
    {
        ROMol ro_mol(*mol);

        std::vector<uint32_t> *invars = new std::vector<uint32_t>(ro_mol.getNumAtoms());
        RDKit::MorganFingerprints::getFeatureInvariants(*mol, *invars);

        auto res = RDKit::MorganFingerprints::getFingerprint(
            *mol, 2,
            invars, nullptr, false,
            true, true, false, nullptr, false);

        auto map = res->getNonzeroElements();

        std::vector<uint64_t> bytes;

        // printf("%d\n", map.size());

        for (auto const &val : map)
        {
            bytes.push_back(val.first);
        }
        delete invars;

        std::vector<uint64_t> *bytes_heap = new std::vector<uint64_t>(bytes);

        return std::unique_ptr<std::vector<uint64_t>>(bytes_heap);
    }

    std::shared_ptr<ExplicitBitVect> copy_explicit_bit_vect(std::shared_ptr<ExplicitBitVect> orig)
    {
        std::shared_ptr<ExplicitBitVect> fingerprint(new ExplicitBitVect(*orig));
        return fingerprint;
    }

    std::unique_ptr<std::vector<uint64_t>> explicit_bit_vect_to_u64_vec(std::shared_ptr<ExplicitBitVect> bitvect)
    {
        std::vector<uint64_t> bytes;
        bytes.reserve(bitvect->dp_bits->num_blocks());
        boost::to_block_range(*bitvect->dp_bits, (std::back_inserter(bytes)));
        std::vector<uint64_t> *bytes_heap = new std::vector<uint64_t>(bytes);
        return std::unique_ptr<std::vector<uint64_t>>(bytes_heap);
    }
}