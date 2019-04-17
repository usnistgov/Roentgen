package gov.nist.microanalysis.roentgen.DataStore;

import java.text.ParseException;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Map.Entry;


import java.util.Set;

import org.apache.commons.math3.geometry.euclidean.oned.Interval;
import org.apache.commons.math3.geometry.partitioning.Region.Location;
import org.apache.commons.math3.util.Pair;

import com.duckandcover.lazy.SimplyLazy;

import gov.nist.juncertainty.UncertainValue;
import gov.nist.microanalysis.roentgen.ArgumentException;
import gov.nist.microanalysis.roentgen.physics.Element;
import gov.nist.microanalysis.roentgen.physics.composition.Composition;
import gov.nist.microanalysis.roentgen.physics.composition.MaterialLabel;
import gov.nist.microanalysis.roentgen.physics.composition.MaterialLabel.MassFraction;

/**
 * @author nicho
 *
 */
public class CompositionFactory implements ICompositionArchive {

	private final Map<String, LazyComposition> mCompositions = new HashMap<>();

	protected abstract static class LazyComposition //
			extends SimplyLazy<Composition> {
		private final String mName;

		private LazyComposition(
				String name
		) {
			mName = name;
		}

	}

	protected static class LazyCompound extends LazyComposition {

		private final String mDefinition;

		LazyCompound(
				String name, String definition
		) {
			super(name);
			mDefinition = definition;
		}

		@Override
		protected Composition initialize() {
			try {
				return Composition.parse(mDefinition);
			} catch (Exception e) {
				e.printStackTrace();
			}
			return null;
		}

	}

	protected static class NotLazy extends LazyComposition {

		private final Composition mComposition;

		NotLazy(
				String name, Composition comp
		) {
			super(name);
			mComposition = comp;
		}

		@Override
		protected Composition initialize() {
			return mComposition;
		}

	}

	private static final LazyComposition K240 = new LazyComposition("K240") {

		@Override
		protected Composition initialize() {
			final Map<Element, Number> res = new HashMap<Element, Number>();
			res.put(Element.Magnesium, new UncertainValue(0.030154, 0.00030154));
			res.put(Element.Silicon, new UncertainValue(0.186986, 0.00186986));
			res.put(Element.Titanium, new UncertainValue(0.059950, 0.00059950));
			res.put(Element.Zinc, new UncertainValue(0.040168, 0.00040168));
			res.put(Element.Zirconium, new UncertainValue(0.074030, 0.00074030));
			res.put(Element.Barium, new UncertainValue(0.268689, 0.000268689));
			res.put(Element.Oxygen, new UncertainValue(0.340023, 0.00340023));
			try {
				return Composition.massFraction("K240", res);
			} catch (ArgumentException e) {
				e.printStackTrace();
			}
			return null;
		}
	};

	/**
	 * K411 EPMA values (SP 260-74)
	 */
	private static final LazyComposition K411e = new LazyComposition("K411e") {

		@Override
		protected Composition initialize() {
			try {
				return Composition.combine("K411e", false, //
						Pair.create(Composition.parse("SiO2"), new UncertainValue(0.5389, 0.0096)), //
						Pair.create(Composition.parse("FeO"), new UncertainValue(0.1448, 0.0027)), //
						Pair.create(Composition.parse("MgO"), new UncertainValue(0.1512, 0.0020)), //
						Pair.create(Composition.parse("CaO"), new UncertainValue(0.1549, 0.0015)));
			} catch (ArgumentException | ParseException e) {
				e.printStackTrace();
			}
			return null;
		}
	};

	/**
	 * K411 certified values (SRM 470)
	 */
	private static final LazyComposition K411 = new LazyComposition("K411") {

		@Override
		protected Composition initialize() {
			try {
				return Composition.combine("K411", false, //
						Pair.create(Composition.parse("SiO2"), new UncertainValue(0.5430, 0.0020)), //
						Pair.create(Composition.parse("FeO"), new UncertainValue(0.1442, 0.0020)), //
						Pair.create(Composition.parse("MgO"), new UncertainValue(0.1467, 0.0020)), //
						Pair.create(Composition.parse("CaO"), new UncertainValue(0.1547, 0.0020)));
			} catch (ArgumentException | ParseException e) {
				e.printStackTrace();
			}
			return null;
		}
	};

	/**
	 * K412 EPMA values (SP 260-74)
	 */
	private static final LazyComposition K412e = new LazyComposition("K412e") {

		@Override
		protected Composition initialize() {
			try {
				return Composition.combine("K412e", false, //
						Pair.create(Composition.parse("SiO2"), new UncertainValue(0.4541, 0.0077)), //
						Pair.create(Composition.parse("FeO"), new UncertainValue(0.0994, 0.0018)), //
						Pair.create(Composition.parse("MgO"), new UncertainValue(0.1966, 0.0025)), //
						Pair.create(Composition.parse("CaO"), new UncertainValue(0.1544, 0.0015)), //
						Pair.create(Composition.parse("Al2O3"), new UncertainValue(0.0934, 0.0029)));
			} catch (ArgumentException | ParseException e) {
				e.printStackTrace();
			}
			return null;
		}
	};

	/**
	 * K412 certified values (SRM 470)
	 */
	private static final LazyComposition K412 = new LazyComposition("K412") {

		@Override
		protected Composition initialize() {
			try {
				return Composition.combine("K412", false, //
						Pair.create(Composition.parse("SiO2"), new UncertainValue(0.4535, 0.0020)), //
						Pair.create(Composition.parse("FeO"), new UncertainValue(0.1010, 0.0020)), //
						Pair.create(Composition.parse("MgO"), new UncertainValue(0.1933, 0.0020)), //
						Pair.create(Composition.parse("CaO"), new UncertainValue(0.1525, 0.0020)), //
						Pair.create(Composition.parse("Al2O3"), new UncertainValue(0.0927, 0.0020)));
			} catch (ArgumentException | ParseException e) {
				e.printStackTrace();
			}
			return null;
		}
	};

	private static final LazyComposition BENITOITE = new LazyCompound("Benitoite", "BaTiSi3O9");
	private static final LazyComposition CAF2 = new LazyCompound("Calcium Fluoride", "CaF2");
	private static final LazyComposition ALUMINA = new LazyCompound("Alumina", "Al2O3");
	private static final LazyComposition SANBORNITE = new LazyCompound("Sanbornite", "BaSi2O5");
	private static final LazyComposition FORSTERITE = new LazyCompound("Forsterite", "Mg2SiO4");
	private static final LazyComposition ZIRCON = new LazyCompound("Zircon", "ZrSiO4");
	private static final LazyComposition ALBITE = new LazyCompound("Albite", "NaAlSi3O8");

	// Must be after all the other static initializers!!!
	private static final CompositionFactory INSTANCE = new CompositionFactory();

	public static CompositionFactory instance() {
		return INSTANCE;
	}

	/**
	 * 
	 */
	private CompositionFactory() {
		mCompositions.put(K240.mName, K240);
		mCompositions.put(BENITOITE.mName, BENITOITE);
		mCompositions.put(CAF2.mName, CAF2);
		mCompositions.put(ALUMINA.mName, ALUMINA);
		mCompositions.put(SANBORNITE.mName, SANBORNITE);
		mCompositions.put(FORSTERITE.mName, FORSTERITE);
		mCompositions.put(K412.mName, K412);
		mCompositions.put(K411.mName, K411);
		mCompositions.put(K412e.mName, K412e);
		mCompositions.put(K411e.mName, K411e);
		mCompositions.put(ZIRCON.mName, ZIRCON);
		mCompositions.put(ALBITE.mName, ALBITE);
	}

	/**
	 *  
	 * @see
	 * gov.nist.microanalysis.roentgen.DataStore.ICompositionArchive#findComposition
	 * (java.lang.String)
	 */
	@Override
	public Curated<Composition> findComposition(
			String name
	) throws Exception {
		final LazyComposition res = mCompositions.getOrDefault(name, null);
		return res != null ? new Curated<Composition>("Factory[" + name + "]", res.get()) : null;
	}

	/**
	 * 
	 * @see gov.nist.microanalysis.roentgen.DataStore.ICompositionArchive#
	 * findMatchingComposition(java.util.Map)
	 */
	@Override
	public Set<Curated<Composition>> findMatchingComposition(
			Map<Element, Interval> criteria
	) throws Exception {
		Set<Curated<Composition>> res = new HashSet<>();
		for(Entry<String, LazyComposition> mesc : mCompositions.entrySet()) {
			Composition comp = mesc.getValue().get();
			boolean matches=true;
			for(Map.Entry<Element, Interval> meei : criteria.entrySet()) {
				final MassFraction mft = MaterialLabel.buildMassFractionTag(comp.getMaterial(), meei.getKey());
				if(meei.getValue().checkPoint(comp.getEntry(mft), 0.0)!=Location.INSIDE) {
					matches=false;
					break;
				}
			}
			if(matches)
				res.add(new Curated<Composition>("Factory[" + comp.toString() + "]", comp));
		}
		return res;
	}

	/**
	 * 
	 * @see gov.nist.microanalysis.roentgen.DataStore.ICompositionArchive#
	 * getCompositionNames()
	 */
	@Override
	public Set<String> getCompositionNames() throws Exception {
		return Collections.unmodifiableSet(mCompositions.keySet());
	}

	/**
	 * 
	 * @see
	 * gov.nist.microanalysis.roentgen.DataStore.ICompositionArchive#add(java.lang.
	 * String, gov.nist.microanalysis.roentgen.physics.composition.Composition)
	 */
	@Override
	public Curated<Composition> add(
			String name, Composition comp
	) {
		final NotLazy value = new NotLazy(comp.toString(), comp);
		mCompositions.put(name, value);
		return new Curated<Composition>("Factory[" + comp.toString() + "]", value.get());
	}

}
