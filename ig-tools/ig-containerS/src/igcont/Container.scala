package igcont

import igcont.trie.Trie
import igcont.kmer.Counter
import igcont.anno.{Record, Anno}
import scala.collection.immutable.HashMap
import scala.collection.mutable

/**
 * Created with IntelliJ IDEA.
 * User: mactep
 * Date: 18.10.13
 * Time: 11:48
 */
class Container(alphabet : String, special : Char, anno_types : Array[String], k : Int) {
  private val _trie  = new Trie()
  private val _kstat = new Counter(alphabet, special, k)
  private val _anno  = new Anno(anno_types)

  def this(alphabet : String, special : Char)= {
    this(alphabet, special, Array.empty, 7)
  }

  def push(seq : String, name : String) : Int = {
    val record = _anno.createRecord(name, seq.size)
    val handle = record.handle()
    val nodes = new mutable.ArrayBuffer[Int](seq.size)
    var key = 0

    // Add to trie
    seq.foreach(c => {
      key = _trie.insert(key, c)
      _trie.setDataOf(key, mutable.ArrayBuffer.empty[(Int, Int)])
      nodes += key
    })

    // Set nodes for annotations
    (0 until nodes.size).foreach(i => {
      record.setNode(i, nodes(i))
      val data =_trie.dataOf(nodes(i)).asInstanceOf[mutable.ArrayBuffer[(Int, Int)]]
      data += ((handle, i))
    })

    // Update kmer index
    _kstat.add(handle, seq, nodes)

    handle
  }

  def record(handle : Int) : Record = _anno.getRecord(handle)

  def record(name : String) : Record = _anno.getRecord(name)

  def seq(handle : Int) : String = {
    seq(record(handle))
  }

  def seq(name : String) : String = {
    seq(record(name))
  }

  def seq(rec : Record) : String = {
    val s = new StringBuilder()
    (0 until rec.size()).foreach(i => {
      s.append(_trie.symbolOf(rec.nodeOf(i)))
    })

    s.toString()
  }

  def data(handle : Int) : Iterable[(Char, HashMap[String, String])] = {
    data(record(handle))
  }

  def data(name : String) : Iterable[(Char, HashMap[String, String])] = {
    data(record(name))
  }

  def data(rec : Record) : Iterable[(Char, HashMap[String, String])] = {
    val result = mutable.ArrayBuffer.fill[(Char, HashMap[String, String])](rec.size())(null)
    (0 until rec.size()).foreach(i => {
      result(i) = (_trie.symbolOf(rec.nodeOf(i)), rec.annotationOf(i))
    })
    result
  }

  def labels() : Iterable[String] = _anno.keys()

  def size() : Int = _anno.size()

  def nodes() : Int = _trie.size()

  // Algorithms

  def find(pattern : String) : Iterable[(String, Int)] = {
    val len = pattern.size - _kstat.k() + 1
    val tmp_trie = new Trie()
    tmp_trie.copyOf(_trie)
    (0 until tmp_trie.size()).foreach(i => tmp_trie.setDataOf(i, false))

    // Get all kmers of the string
    val kmers = (0 until len).map(i => pattern slice(i, i + _kstat.k()))

    // Nodes of last kmer
    var last : Iterable[Int] = null
    // Handles of records in result
    var handles : mutable.Set[Int] = null

    kmers.foreach(kmer => {
      val g = _kstat.get(kmer)
      g match {
        case Some((lst, set)) => {
          last = lst
          lst.foreach(n => tmp_trie.setDataOf(n, true))
          if (handles == null) {
            handles = set
          }
          else {
            handles = handles intersect set
          }
        }
        case None => return Array.empty[(String, Int)]
      }
    })

    val result = mutable.ArrayBuffer.empty[(String, Int)]
    last.foreach(node => {
      var tmp = node
      var counter = 0
      var colored = true
      while (counter != len && colored) {
        counter += 1
        tmp = tmp_trie.parentOf(tmp)
        colored = tmp_trie.dataOf(tmp).asInstanceOf[Boolean]
      }

      if (counter == len) {
        _trie.dataOf(node).asInstanceOf[mutable.ArrayBuffer[(Int, Int)]].foreach(pair => {
          // Special guard to choose only right sequences
          // Filter is the set of sequences, each of that has all the kmers of pattern
          if (handles.contains(pair._1)) {
            result += ((_anno.getRecord(pair._1).name(), pair._2 - len + 1))
          }
        })
      }
    })

    result
  }
}
