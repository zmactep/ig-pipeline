package igcont.anno

import scala.collection.immutable
import scala.collection.mutable

/**
 * Created with IntelliJ IDEA.
 * User: mactep
 * Date: 17.10.13
 * Time: 17:36
 */
class Anno(anno_types : Array[String]) {
  private val _cont   = new mutable.ArrayBuffer[Record]()
  private val _names  = new mutable.HashMap[String, Record]()
  private val _atypes = immutable.HashMap(anno_types.map(i => (i, mutable.ArrayBuffer.empty[String])).toSeq:_*)


  def createRecord(name : String, size : Int) : Record = addRecord(new Record(_cont.size, name, size, _atypes))

  def getRecord(i : Int) : Record = _cont(i)

  def getRecord(name : String) : Record = _names.get(name).get

  def annotationOf(rec : Int, pos : Int) : immutable.HashMap[String, String] = _cont(rec).annotationOf(pos)

  def nodeOf(rec : Int, pos : Int) : Int = _cont(rec).nodeOf(pos)

  def check(atype : String, aval : String) : Boolean = _atypes(atype).contains(aval)

  def add(atype : String, aval : String) = _atypes(atype) += aval

  def keys : Iterable[String] = _atypes.keys

  def valuesOf(atype : String) : Iterable[String] = _atypes(atype)

  def size : Int = _cont.size

  def fullsize : Int = _cont.foldRight(0)((rec, acc) => acc + rec.size)

  def types : immutable.HashMap[String, mutable.ArrayBuffer[String]] = _atypes

  private def addRecord(rec : Record) : Record = {
    _cont += rec
    _names.put(rec.name, rec)
    rec
  }
}
