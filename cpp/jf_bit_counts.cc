// Take in a Jellyfish hash, and replace the count of all kmers with 1.
// This essentially reduces a hash to binary values indicating whether or not a kmer
// is present.
// Cobbled together from various Jellyfish subcommand code.

#include <iostream>
#include <fstream>


#include <jellyfish/jellyfish.hpp>


namespace err = jellyfish::err;

template<typename iterator, typename writer_type>
void bit_count(iterator& it, std::ostream &out, writer_type& writer) {
    while (it.next()) {
        writer.write(out, it.key(), 1);
    }
}

int main(int argc, char *argv[])
{
    if (argc != 3)
        err::die(err::msg() << "Usage: " << argv[0] << " in.jf out.jf");

    std::ifstream is(argv[1]);
    if (!is.good()) {
        err::die(err::msg() << "Could not read input file: " << argv[1]);
    }

    std::ofstream out(argv[2]);
    if (!out.good()) {
        err::die(err::msg() << "Could not open output file: " << argv[2]);
    }

    jellyfish::file_header header;
    header.read(is);
    jellyfish::mer_dna::k(header.key_len() / 2);

    if (!header.format().compare(binary_dumper::format)) {
        binary_reader reader(is, &header);
        header.write(out);
        binary_writer writer(header.counter_len(), header.key_len());
        bit_count(reader, out, writer);
    } else if (!header.format().compare(text_dumper::format)) {
        err::die(err::msg() << "bit_counts does not support Jellyfish text format");
    }

    out.close();

    return 0;
}