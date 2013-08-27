#include <vector>
#include <log4cxx/logger.h>

#include "database.h"
#include "generic_exception.h"
#include "reverse_complement.h"
#include "kmer_generator.h"
#include "fasta_reader.h"
#include "read_processor.h"

extern log4cxx::LoggerPtr logger;

void DatabaseFiller::insert2db(std::string * kmer, std::string * sequence){

#pragma omp critical
	{
		if (kmer2listOfSeq->end() == kmer2listOfSeq->find(kmer)) {
			std::set<std::string *, Compare> source;
			source.insert(sequence);
			kmer2listOfSeq->insert(std::make_pair(kmer, source));
		} else {
			(*kmer2listOfSeq)[kmer].insert(sequence);
		}
	}
}

bool DatabaseFiller::operator()(const Read &r) {
	//TODO(Feodorov) do we need reverse complement reads and k-mers?
	std::string * name = new std::string(r.getName());
	std::string * sequence = new std::string(r.getSeq());

#pragma omp critical
	{
		name2seq->insert(std::make_pair(name, sequence));
		seq2name->insert(std::make_pair(sequence, name));
	}

	std::vector<std::string> kmers = KmerGenerator::getKmers(r, kmer_size);
	for (std::vector<std::string>::const_iterator it = kmers.begin(); it != kmers.end(); ++it) {
		insert2db(new std::string(*it), sequence);
	}
	return false;
}

Database::Database(const struct settings_t& settings) {
	FastaReader reader(settings.reference_file);
	DatabaseFiller filler(settings.kmer_size);
	ReadProcessor rp(settings.max_threads);
	rp.readAndProcessSingleThread(reader, filler);

	LOG4CXX_DEBUG(logger, "Number of reference reads is " << rp.get_num_processed_reads());
	name2seq = filler.getName2seq();
	seq2name = filler.getSeq2name();
	kmer2listOfSeq = filler.getKmer2listOfSeq();
}

Database::~Database() {
	for (auto it = name2seq->begin(); it != name2seq->end(); ++it) {
		delete it->first;
		delete it->second;
	}

	for (auto it = kmer2listOfSeq->begin(); it != kmer2listOfSeq->end(); ++it) {
		delete it->first;
	}

	delete name2seq;
	delete seq2name;
	delete kmer2listOfSeq;
}

void Database::get_sequence_by_name(const std::string& name, std::string& out_seq) const {
	auto it = name2seq->find(const_cast<std::string *>(&name));
	if (name2seq->end() == it) {
		throw GenericException("Element not found: " + name);
	}
	out_seq.assign(*(it->second));
}

void Database::get_name_by_sequence(const std::string& seq, std::string& out_name) const {
	auto it = seq2name->find(const_cast<std::string *>(&seq));
	if (name2seq->end() == it) {
		throw GenericException("Element not found: " + seq);
	}
	out_name.assign(*(it->second));
}

void Database::get_sequences_for_kmer(const std::string& kmer, std::set<std::string *, Compare>& out_seq) const {
	auto it = kmer2listOfSeq->find(const_cast<std::string *>(&kmer));
	if (kmer2listOfSeq->end() != it) {
		out_seq.insert(it->second.begin(), it->second.end());
	}
}

int Database::get_num_sequences() const {
	return name2seq->size();
}

int Database::get_num_kmers() const {
	return kmer2listOfSeq->size();
}

std::map<std::string *, std::string *>::const_iterator Database::get_data_iterator() const {
	return name2seq->begin();
}

std::map<std::string *, std::set<std::string *, Compare>, Compare>::const_iterator Database::get_kmer_iterator() const {
	return kmer2listOfSeq->begin();
}
