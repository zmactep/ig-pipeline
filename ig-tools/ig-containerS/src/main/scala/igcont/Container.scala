package igcont

import igcont.trie.Trie
import igcont.anno.{Record, Anno}
import scala.collection.immutable.HashMap
import scala.collection.mutable
import alicont.{AlicontFactory, AbstractAlicont, AlignmentResult}
import scala.collection.mutable.ArrayBuffer
import igcont.kmer.bit.Counter

import alicont.algorithms.AlgorithmType.AlgorithmType

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
  private var _depth = 0

  def this(alphabet : String, special : Char)= {
    this(alphabet, special, Array.empty, 7)
  }

  def push(seq : String, name : String) : Int = {
    val record = _anno.createRecord(name, seq.size)
    val handle = record.handle
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

    // Update max depth
    if (_depth < seq.size) {
      _depth = seq.size
    }

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
    (0 until rec.size).foreach(i => {
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
    val result = mutable.ArrayBuffer.fill[(Char, HashMap[String, String])](rec.size)(null)
    (0 until rec.size).foreach(i => {
      result(i) = (_trie.symbolOf(rec.nodeOf(i)), rec.annotationOf(i))
    })
    result
  }

  def labels : Iterable[String] = _anno.keys

  def size : Int = _anno.size

  def fullsize : Int = _anno.fullsize

  def nodes : Int = _trie.size

  def addAnnotations(atype : String, annotations : Array[String]) : Unit = {
    annotations.foreach(a => if (!_anno.check(atype, a)) _anno.add(atype, a))
  }

  // Algorithms

  def find(pattern : String) : Iterable[(String, Int)] = {
    val len = pattern.size - _kstat.k + 1
    val tmp_trie = new Trie()
    tmp_trie.copyOf(_trie)
    (0 until tmp_trie.size).foreach(i => tmp_trie.setDataOf(i, false))

    // Get all kmers of the string
    val kmers = (0 until len).map(i => pattern slice(i, i + _kstat.k))

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
        tmp_trie.parentOf(tmp) match {
          case Some(node_id) => {
            tmp = node_id
            colored = tmp_trie.dataOf(tmp).asInstanceOf[Boolean]
            counter += 1
          }
          case None => colored = false
        }
      }

      if (counter == len) {
        _trie.dataOf(node).asInstanceOf[mutable.ArrayBuffer[(Int, Int)]].foreach(pair => {
          // Special guard to choose only right sequences
          // Filter is the set of sequences, each of that has all the kmers of pattern
          if (handles.contains(pair._1)) {
            result += ((_anno.getRecord(pair._1).name, pair._2 - len + 1))
          }
        })
      }
    })

    result
  }

  def alignment_template(alicont : AbstractAlicont, callback : (AlignmentResult) => Unit) : Unit = {
    val target = new mutable.StringBuilder()
    val fork_stack = new mutable.Stack[Int]()
    var from_leaf = false

    fork_stack.push(0)

    var node = _trie.nextOf(0)
    while (node != 0) {
      // Check if we made a DFS-jump
      if (from_leaf) {
        val parent = _trie.parentOf(node).getOrElse(0)
        while (fork_stack.top != parent) {
          fork_stack.pop()
          alicont.pop()
        }
      }

      from_leaf = false
      target.append(_trie.symbolOf(node))
      if (_trie.isFork(node)) {
        fork_stack.push(node)
        alicont.push(target.toString())
        target.clear()
      }
      if (_trie.isLeaf(node)) {
        from_leaf = true
        fork_stack.push(node)

        alicont.push(target.toString())
        target.clear()

        val (score, (q, t)) = alicont.alignment()
        val node_data = _trie.dataOf(node).asInstanceOf[ArrayBuffer[(Int, Int)]]
        node_data.foreach(tpl => {
          val align = new AlignmentResult(score, q, t, (seq(tpl._1), record(tpl._1)))
          // Your logic here
          callback(align)
        })
      }

      node = _trie.nextOf(node)
    }
  }

  def alignment_inner(alicont : AbstractAlicont, n : Int) : Iterable[AlignmentResult] = {
    val result = new mutable.PriorityQueue[AlignmentResult]()(Ordering.by(a => -a.score))

    def n_callback(align : AlignmentResult) : Unit = {
      if (result.size == n) {
        val m = result.head
        if (align.score > m.score) {
          result += align
          result.dequeue()
        }
      }
      else {
        result += align
      }
    }

    alignment_template(alicont, n_callback)

    result.dequeueAll.reverse.toIterable
  }

  def alignment_inner(alicont : AbstractAlicont, prct : Double) : Iterable[AlignmentResult] = {
    val result = new mutable.PriorityQueue[AlignmentResult]()(Ordering.by(a => a.score))

    def n_callback(align : AlignmentResult) : Unit = {
      if (align.similarity >= prct) {
        result += align
      }
    }

    alignment_template(alicont, n_callback)

    result.dequeueAll.reverse.toIterable
  }

  def alignment(query : String, gap : Int, score_matrix : Array[Array[Int]],
                algo_type : AlgorithmType, n : Int) : Iterable[AlignmentResult] =
    alignment_inner(AlicontFactory.createSimpleAlicont(_depth, query, gap, score_matrix, algo_type), n)

  def alignment(query : String, gap : Int, score_matrix : Array[Array[Int]],
                algo_type : AlgorithmType, prct : Double) : Iterable[AlignmentResult] =
    alignment_inner(AlicontFactory.createSimpleAlicont(_depth, query, gap, score_matrix, algo_type), prct)

  def affine_alignment(query : String, gap_open : Int, gap_ext : Int,
                       score_matrix : Array[Array[Int]], algo_type : AlgorithmType, n : Int) : Iterable[AlignmentResult] =
    alignment_inner(AlicontFactory.createAffineAlicont(_depth, query, gap_open, gap_ext, score_matrix, algo_type), n)

  def affine_alignment(query : String, gap_open : Int, gap_ext : Int,
                       score_matrix : Array[Array[Int]], algo_type : AlgorithmType, prct : Double) : Iterable[AlignmentResult] =
    alignment_inner(AlicontFactory.createAffineAlicont(_depth, query, gap_open, gap_ext, score_matrix, algo_type), prct)


  def annotate(query : String, gap : Int, score_matrix : Array[Array[Int]],
               algo_type : AlgorithmType, n : Int) : Iterable[(Char, HashMap[String, String])] = {
    val result = ArrayBuffer.fill[(Char, HashMap[String, String])](query.size)(null)
    val align_result = alignment(query, gap, score_matrix, algo_type, n)
    val preresult = ArrayBuffer.fill[ArrayBuffer[HashMap[String, String]]](query.size)(ArrayBuffer.empty)

    def annotate_by_one(anno : AlignmentResult) = {
      var i = 0
      anno.get.zipWithIndex.foreach(tpl => {
        val ((_, a), j) = tpl
        if (anno.query(j) != '-') {
          if (a != null) {
            preresult(i) += a
          }
          i += 1
        }
      })
    }

    // Get preresult
    align_result.foreach(align => {
      annotate_by_one(align)
    })

    // Merge
    val keys = _anno.keys
    val tmp_hash = mutable.HashMap.empty[String, Int]
    val anno_builder = mutable.HashMap.empty[String, String]

    // For each letter
    preresult.zipWithIndex.foreach(tpl => {
      val (lst, i) = tpl

      anno_builder.clear()
      // For each annotation type
      keys.foreach(key => {
        var max_elem : String = null
        var max_val  = 0
        tmp_hash.clear()
        // For each annotation value in list
        lst.foreach(elem => {
          if (elem.contains(key)) {
            val eval = elem(key)
            val count = tmp_hash.getOrElse(eval, 0) + 1
            tmp_hash(elem(key)) = count
            if (count > max_val && !eval.isEmpty) {
              max_val = count
              max_elem = eval
            }
          }
        })

        // We have an annotation value for key
        if (max_elem != null) {
          anno_builder(key) = max_elem
        }
      })

      result(i) = (query(i), HashMap[String, String](anno_builder.toSeq:_*))
    })

    result
  }
}
