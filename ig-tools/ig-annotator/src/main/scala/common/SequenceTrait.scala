package common

/**
 * Created with IntelliJ IDEA.
 * User: mactep
 * Date: 31.10.13
 * Time: 12:11
 */
class SequenceTrait(t : SequenceType.SequenceType) {
  private val _type = t


  def alphabet : String = {
    _type match {
      case SequenceType.NUCLEO => "ACGTN"
      case SequenceType.AMINO  => "ARNDCEQGHILKMFPSTWYVX"
    }
  }

  def special : Char = {
    _type match {
      case SequenceType.NUCLEO => 'N'
      case SequenceType.AMINO  => 'X'
    }
  }

  def k : Int = {
    _type match {
      case SequenceType.NUCLEO => 7
      case SequenceType.AMINO  => 3
    }
  }
}
