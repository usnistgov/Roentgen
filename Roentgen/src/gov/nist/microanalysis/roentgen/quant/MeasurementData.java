package gov.nist.microanalysis.roentgen.quant;

import java.text.DecimalFormat;

import org.apache.commons.math3.geometry.euclidean.threed.Vector3D;

import com.duckandcover.html.IToHTML;
import com.duckandcover.html.Table;
import com.google.common.base.Objects;
import com.google.common.base.Optional;

import gov.nist.microanalysis.roentgen.math.uncertainty.UncertainValue;
import gov.nist.microanalysis.roentgen.physics.composition.Composition;
import gov.nist.microanalysis.roentgen.physics.composition.Layer;

public class MeasurementData
   implements
   IToHTML {

   private final UncertainValue mLiveTime;
   private final UncertainValue mProbeCurrent;
   private final UncertainValue mBeamEnergy;
   private Vector3D mSurfaceOrientation;
   private Optional<Layer> mCoating;
   private final Optional<Composition> mComposition;

   public final static Vector3D BEAM_NORMAL_TO_SURFACE = Vector3D.MINUS_K;

   private MeasurementData(final UncertainValue liveTime, final UncertainValue probeCurrent, final UncertainValue beamEnergy, final Composition comp) {
      assert (liveTime != null) && (liveTime.doubleValue() > 0.0);
      assert (probeCurrent != null) && (probeCurrent.doubleValue() > 0.0);
      assert (beamEnergy != null) && (beamEnergy.doubleValue() >= 0.1);
      mLiveTime = liveTime;
      mProbeCurrent = probeCurrent;
      mBeamEnergy = beamEnergy;
      mComposition = Optional.fromNullable(comp);
      mCoating = Optional.absent();
      mSurfaceOrientation = BEAM_NORMAL_TO_SURFACE;
   }

   @Override
   public int hashCode() {
      return Objects.hashCode(mBeamEnergy, mCoating, mComposition, mLiveTime, mProbeCurrent, mSurfaceOrientation);
   }

   @Override
   public boolean equals(final Object obj) {
      if(this == obj)
         return true;
      if(obj == null)
         return false;
      if(getClass() != obj.getClass())
         return false;
      final MeasurementData other = (MeasurementData) obj;
      if(!mBeamEnergy.equals(other.mBeamEnergy))
         return false;
      if(!mCoating.equals(other.mCoating))
         return false;
      if(!mComposition.equals(other.mComposition))
         return false;
      if(!mLiveTime.equals(other.mLiveTime))
         return false;
      if(!mProbeCurrent.equals(other.mProbeCurrent))
         return false;
      if(!mSurfaceOrientation.equals(other.mSurfaceOrientation))
         return false;
      return true;
   }

   /**
    * Creates a {@link MeasurementData} object that is suitable for use as an
    * unknown.
    *
    * @param liveTime
    * @param probeCurrent
    * @param beamEnergy
    * @return MeasurementDatum
    */
   static public MeasurementData unknown(final UncertainValue liveTime, final UncertainValue probeCurrent, final UncertainValue beamEnergy) {
      final MeasurementData res = new MeasurementData(liveTime, probeCurrent, beamEnergy, null);
      assert res.isSuitableAsUnknown();
      return res;
   }

   /**
    * Creates a {@link MeasurementData} object that is suitable for use as a
    * standard.
    *
    * @param liveTime
    * @param probeCurrent
    * @param beamEnergy
    * @param comp
    * @return MeasurementDatum
    */
   static public MeasurementData standard(final UncertainValue liveTime, final UncertainValue probeCurrent, final UncertainValue beamEnergy, final Composition comp) {
      final MeasurementData res = new MeasurementData(liveTime, probeCurrent, beamEnergy, null);
      assert res.isSuitableAsUnknown();
      return res;
   }

   /**
    * Modifies a {@link MeasurementData} object to specify a coating.
    *
    * @param coating
    * @return {@link MeasurementData}
    */
   public MeasurementData withCoating(final Layer coating) {
      mCoating = Optional.fromNullable(coating);
      return this;
   }

   /**
    * Modifies a {@link MeasurementData} object to specify an orientation of the
    * sample relative to the beam.
    *
    * @param orient
    * @return {@link MeasurementData}
    */
   public MeasurementData withOrientation(final Vector3D normal) {
      mSurfaceOrientation = normal.normalize();
      return this;
   }

   /**
    * Modifies a {@link MeasurementData} object to specify that the sample is
    * oriented normal to the beam.
    *
    * @return {@link MeasurementData}
    */
   public MeasurementData normalToBeam() {
      mSurfaceOrientation = BEAM_NORMAL_TO_SURFACE;
      return this;
   }

   /**
    * Does the object contain all the data necessary to be used as an unknown?
    *
    * @return boolean
    */
   public boolean isSuitableAsUnknown() {
      assert (mLiveTime != null) && (mLiveTime.doubleValue() > 0.0);
      assert (mProbeCurrent != null) && (mProbeCurrent.doubleValue() > 0.0);
      assert (mBeamEnergy != null) && (mBeamEnergy.doubleValue() >= 0.1);
      return true;
   }

   /**
    * Does the object contain all the data necessary to be used as a standard?
    *
    * @return
    */
   public boolean isSuitableAsStandard() {
      return isSuitableAsUnknown() && mComposition.isPresent();
   }

   /**
    * Is the sample coated?
    *
    * @return boolean
    */
   public boolean isCoated() {
      return mCoating.isPresent();
   }

   /**
    * @return An {@link UncertainValue} in seconds
    */
   public UncertainValue getLiveTime() {
      return mLiveTime;
   }

   /**
    * @return An {@link UncertainValue} in nA
    */
   public UncertainValue getProbeCurrent() {
      return mProbeCurrent;
   }

   /**
    * @return An {@link UncertainValue} in keV
    */
   public UncertainValue getBeamEnergy() {
      return mBeamEnergy;
   }

   /**
    * Returns the optional composition of the material from which the
    * measurement was taken. Not typically present for an unknown but must be
    * present for a standard.
    *
    * @return Optional&lt;Composition&gt;
    */
   public Optional<Composition> getComposition() {
      return mComposition;
   }

   /**
    * Returns the orientation of the sample relative to the beam which is
    * assumed to propagate down towards larger Z.
    *
    * @return Vector3D (normalized to unity length)
    */
   public Vector3D getSurfaceNormal() {
      return mSurfaceOrientation;
   }

   /**
    * The sample is oriented such that the electron beam is normal to the
    * surface.
    *
    * @return boolean
    */
   public boolean isNormal() {
      return mSurfaceOrientation.equals(BEAM_NORMAL_TO_SURFACE);
   }

   public Optional<Layer> getCoating() {
      return mCoating;
   }

   @Override
   public String toHTML(final Mode mode) {
      if(mode == Mode.TERSE) {
         return mComposition.isPresent() ? mComposition.get().toHTML(mode) : "Unknown";
      } else {
         final Table res = new Table();
         final DecimalFormat df = new DecimalFormat("0.0");
         res.addRow(Table.th("Item"), Table.th("Units"), Table.th("Value"));
         res.addRow(Table.td("Live time"), Table.td("(s)"), Table.td(mLiveTime.toHTML(mode)));
         res.addRow(Table.td("Probe current"), Table.td("(nA)"), Table.td(mProbeCurrent.toHTML(mode)));
         res.addRow(Table.td("Beam energy"), Table.td("(keV)"), Table.td(mProbeCurrent.toHTML(mode)));
         if(isNormal())
            res.addRow(Table.td("Surface angle"), Table.td("-"), Table.td("Normal to beam"));
         else
            res.addRow(Table.td("Surface normal"), Table.td("-"), Table.td(mSurfaceOrientation.toString(df)));
         if(mCoating.isPresent())
            res.addRow(Table.td("Surface coating"), Table.td("-"), Table.td(mCoating.get().toHTML(mode)));
         else
            res.addRow(Table.td("Surface coating"), Table.td("-"), Table.td("None"));
         if(mComposition.isPresent())
            res.addRow(Table.td("Composition"), Table.td("-"), Table.td(mComposition.get().toHTML(mode)));
         return res.toHTML(mode);
      }
   }
}
