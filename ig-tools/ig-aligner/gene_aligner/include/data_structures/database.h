#ifndef DATABASE_H_
#define DATABASE_H_

#include <string>
#include <map>
#include <set>
#include "comparator.h"
#include "read.h"
#include "settings.h"

class DatabaseFiller {
public:
	DatabaseFiller(int kmer_size) : kmer_size(kmer_size),
						name2seq(new std::map<std::string *, std::string *, Compare>()),
						seq2name(new std::map<std::string *, std::string *, Compare>()),
						kmer2listOfSeq(new std::map<std::string *, std::set<std::string *, Compare>, Compare>) {};
	bool operator()(const Read &r);

	std::map<std::string *, std::string *, Compare> * getName2seq() const {
		return name2seq;
	}
	std::map<std::string *, std::string *, Compare> * getSeq2name() const {
		return seq2name;
	}
	std::map<std::string *, std::set<std::string *, Compare>, Compare> * getKmer2listOfSeq() const {
		return kmer2listOfSeq;
	}

private:
	void insert2db(std::string * kmer, std::string * source_sequence);
	std::map<std::string *, std::string *, Compare> * name2seq;
	std::map<std::string *, std::string *, Compare> * seq2name;
	std::map<std::string *, std::set<std::string *, Compare>, Compare> * kmer2listOfSeq;
	const int kmer_size;
};

class Database {
public:
	Database(const struct settings_t& settings);
	~Database();
	int get_num_sequences() const;
	int get_num_kmers() const;
	void get_sequence_by_name(const std::string& name, std::string & out_seq) const;
	void get_name_by_sequence(const std::string& seq, std::string & out_name) const;
	void get_sequences_for_kmer(const std::string& kmer, std::set<std::string *, Compare> & out_seq) const;
	std::map<std::string *, std::string *>::const_iterator get_data_iterator() const;
	std::map<std::string *, std::set<std::string *, Compare>, Compare>::const_iterator get_kmer_iterator() const;
private:
	std::map<std::string *, std::string *, Compare> * name2seq;
	std::map<std::string *, std::string *, Compare> * seq2name;
	std::map<std::string *, std::set<std::string *, Compare>, Compare> * kmer2listOfSeq;
};


#endif /* DATABASE_H_ */
