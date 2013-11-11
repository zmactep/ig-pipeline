package igcont.kmer.slow

import scala.collection.immutable.TreeMap
import scala.collection.mutable

/**
 * Created with IntelliJ IDEA.
 * User: mactep
 * Date: 07.10.13
 * Time: 17:03
 */
class Counter(a : String, s : Char, ksize : Int) {
  private val _kmers   : mutable.HashMap[Int, mutable.HashSet[Int]] = new mutable.HashMap[Int, mutable.HashSet[Int]]
  private val _labels  : mutable.HashMap[Int, mutable.HashSet[Int]] = new mutable.HashMap[Int, mutable.HashSet[Int]]
  private val _alpha   : TreeMap[Char, Int] = TreeMap[Char, Int](a.zipWithIndex.toSeq:_*)
  private val _special : Char               = s
  private val _ksize   : Int                = ksize

  def add(handle : Int, seq : String, i : Iterable[Int]) : Unit = {
    val kmers = hashs(seq)
    (kmers zip i).foreach(tpl => {
      val (k, n) = tpl
      if (!_kmers.contains(k)) {
        _kmers.put(k, new mutable.HashSet[Int]())
        _labels.put(k, new mutable.HashSet[Int]())
      }
      _kmers.get(k).get.add(n)
      _labels.get(k).get.add(handle)
    })
  }

  def check(kmer : String, use_special : Boolean) : Boolean =
    kmer.forall(c => _alpha.contains(c) || (use_special && c == _special))

  def hashs(seq : String) : Iterable[Int] = {
    val a = new Hasher(_ksize, _alpha)
    seq.foreach(c => a.add(c))
    a.get
  }

  def hash(kseq : String) : Int =
    if (check(kseq, use_special = false)) hashs(kseq).head else -1

  def get(kmer : String) : Option[(Iterable[Int], mutable.Set[Int])] = {
    if (check(kmer, use_special = true)) {
      val specs = kmer.count(c => c == _special)

      if (specs == 1) {
        val result_n = new mutable.HashSet[Int]()
        val result_h = new mutable.HashSet[Int]()
        _alpha.foreach(c => {
          val rep = hash(kmer.replace(_special, c._1))
          if (_kmers.contains(rep)) {
            result_n ++= _kmers.get(rep).get
            result_h ++= _labels.get(rep).get
          }
        })
        if (result_n.size > 0) {
          return Some(result_n, result_h)
        }
      }

      if (specs == 0) {
        val rep = hash(kmer)
        if (_kmers.contains(rep)) {
          return Some(_kmers.get(rep).get, _labels.get(rep).get)
        }
      }
    }
    None
  }

  def k : Int = {
    _ksize
  }

  def special : Char = _special
}
