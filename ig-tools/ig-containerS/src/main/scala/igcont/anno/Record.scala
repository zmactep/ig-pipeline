package igcont.anno

import scala.collection.mutable
import scala.collection.immutable

import alicont.common

/**
 * Created with IntelliJ IDEA.
 * User: mactep
 * Date: 17.10.13
 * Time: 17:37
 */
class Record(id : Int, n : String, s : Int,
             anno_types : immutable.HashMap[String, mutable.ArrayBuffer[String]]) extends common.Record {
  private val _id     = id
  private val _name   = n
  private val _cont   = mutable.ArrayBuffer.fill[Annotation](s)(null)
  private val _atypes = anno_types


  private[igcont] def setNode(i : Int, node : Int) : Unit = {
    _cont(i) = new Annotation(node, _atypes.size)
  }

  def setAnnotation(pos : Int, atype : String, aval : String) : (Int, Int) = {
    var iaval = -1
    val iatype = atype2idx(atype).getOrElse(-1)

    if (_cont(pos) == null) {
      println("Warning! Cannot make an annotation without node! Use setNode(i, node) first. Default initialization: 0.")
      setNode(pos, 0)
    }

    aval2idx(atype, aval) match {
      case Some(i : Int) => {
        _cont(pos).set(iatype, i)
        iaval = i
      }
      case None => {
        iaval = _atypes(atype).size
        _atypes(atype) += aval
        _cont(pos).set(iatype, iaval)
      }
    }

    (iatype, iaval)
  }

  def setAnnotation(pos : Int, atype : Int, aval : Int) : (Int, Int) = {
    if (pos != -1) {
      if (_cont(pos) == null) {
        println("Warning! Cannot make an annotation without node! Use setNode(i, node) first. Default initialization: 0.")
        setNode(pos, 0)
      }

      _cont(pos).set(atype, aval)
    }
    (atype, aval)
  }

  def name : String = _name

  def handle : Int = _id

  def size : Int = _cont.size

  def annotationOf(pos : Int) : immutable.HashMap[String, String] = {
    if (_cont(pos) == null) {
      println("Warning! Cannot get an annotation without node! Use setNode(i, node) first.")
      return immutable.HashMap.empty[String, String]
    }
    immutable.HashMap(_atypes.keys.zipWithIndex.map(k =>
      (k._1, idx2aval(k._1, _cont(pos).dataOf(k._2)).getOrElse("")))
      .toSeq:_*)
  }

  private[igcont] def nodeOf(i : Int) : Int = _cont(i).node

  private def atype2idx(s : String) : Option[Int] = {
    if (_atypes.keys.toSeq.contains(s)) {
      return Some(_atypes.keys.toSeq.indexOf(s))
    }
    None
  }

  private def idx2atype(i : Int) : Option[String] = {
    if (_atypes.keys.size < i) {
      Some(_atypes.keys.toSeq(i))
    }
    None
  }

  private def aval2idx(atype: String, aval : String) : Option[Int] = {
    val sz = _atypes(atype).size
    if (sz > 0) {
      (0 until sz).foreach(i =>
        if (_atypes(atype)(i) == aval)
          return Some(i)
      )
    }
    None
  }

  private def idx2aval(atype : String, i : Int) : Option[String] = {
    val vec = _atypes(atype)
    if (i >= 0 && i < vec.size) {
      return Some(vec(i))
    }
    None
  }
}
