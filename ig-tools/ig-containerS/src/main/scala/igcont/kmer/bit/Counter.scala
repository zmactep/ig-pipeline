package igcont.kmer.bit

import scala.collection.mutable
import scala.collection.immutable.TreeMap

/**
 * Created with IntelliJ IDEA.
 * User: mactep
 * Date: 11.11.13
 * Time: 14:08
 */
class Counter(a : String, s : Char, ksize : Int) {
  private val _kmers   : mutable.HashMap[Int, mutable.HashSet[Int]] = new mutable.HashMap[Int, mutable.HashSet[Int]]
  private val _labels  : mutable.HashMap[Int, mutable.HashSet[Int]] = new mutable.HashMap[Int, mutable.HashSet[Int]]
  private val _alpha   : TreeMap[Char, Int] = TreeMap[Char, Int](a.zipWithIndex.toSeq:_*)
  private val _special : Char               = s
  private val _ksize   : Int                = ksize
  private val _hasher  : Hasher          = new Hasher(_ksize, _alpha)

  def add(handle : Int, seq : String, i : Iterable[Int]) : Unit = {
    var k = hash(seq)
    var j = _ksize
    i.foreach(n => {
      if (!_kmers.contains(k)) {
        _kmers.put(k, new mutable.HashSet[Int]())
        _labels.put(k, new mutable.HashSet[Int]())
      }
      _kmers.get(k).get.add(n)
      _labels.get(k).get.add(handle)
      if (j < seq.size) {
        k = _hasher.add(seq(j))
        j += 1
      }
    })
  }

  def check(kmer : String, use_special : Boolean) : Boolean =
    kmer.forall(c => _alpha.contains(c) || (use_special && c == _special))

  def hash(kseq : String) : Int = {
    if (check(kseq, use_special = false)) {
      var i = -1
      _hasher.clear()
      kseq.zipWithIndex.takeWhile(tpl => tpl._2 < _ksize).foreach(tpl => i = _hasher.add(tpl._1))
      return i
    }
    -1
  }

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
