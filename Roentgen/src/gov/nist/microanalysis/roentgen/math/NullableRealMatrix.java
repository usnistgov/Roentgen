package gov.nist.microanalysis.roentgen.math;

import org.apache.commons.math3.exception.NotStrictlyPositiveException;
import org.apache.commons.math3.exception.OutOfRangeException;
import org.apache.commons.math3.linear.AbstractRealMatrix;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;

/**
 * Creates a matrix in which one or both of the dimensions can be zero.
 *
 *
 * @author Nicholas W. M. Ritchie
 *
 */
public class NullableRealMatrix extends AbstractRealMatrix {

	public static RealMatrix build(final int rows, final int cols) {
		if ((rows == 0) || (cols == 0))
			return new NullableRealMatrix(rows, cols);
		else
			return MatrixUtils.createRealMatrix(rows, cols);
	}

	public static RealMatrix build(final RealMatrix rm) {
		if ((rm.getRowDimension() == 0) || (rm.getColumnDimension() == 0))
			return new NullableRealMatrix(rm);
		else
			return rm.copy();
	}

	public static RealMatrix build(final double[][] vals) {
		if ((vals.length == 0) || (vals[0].length == 0))
			return new NullableRealMatrix(vals);
		else
			return MatrixUtils.createRealMatrix(vals);
	}

	final int mRowDim;
	final int mColDim;

	private NullableRealMatrix(final int rows, final int cols) {
		assert (rows == 0) || (cols == 0);
		mRowDim = Math.max(0, rows);
		mColDim = Math.max(0, cols);
	}

	private NullableRealMatrix(final RealMatrix rm) {
		mRowDim = rm.getRowDimension();
		mColDim = rm.getColumnDimension();
	}

	private NullableRealMatrix(final double[][] vals) {
		mRowDim = vals.length;
		mColDim = vals.length > 0 ? vals[0].length : 0;
	}

	private NullableRealMatrix(final int rows, final int cols, final RealMatrix matrix) {
		mRowDim = rows;
		mColDim = cols;
	}

	@Override
	public int getRowDimension() {
		return mRowDim;
	}

	@Override
	public int getColumnDimension() {
		return mColDim;
	}

	@Override
	public RealMatrix createMatrix(final int rowDimension, final int columnDimension) //
			throws NotStrictlyPositiveException {
		return new NullableRealMatrix(rowDimension, columnDimension);
	}

	@Override
	public RealMatrix copy() {
		return new NullableRealMatrix(mRowDim, mColDim);
	}

	@Override
	public double getEntry(final int row, final int column) throws OutOfRangeException {
		throw new OutOfRangeException(0, 0, -1);
	}

	@Override
	public void setEntry(final int row, final int column, final double value) throws OutOfRangeException {
		throw new OutOfRangeException(0, 0, -1);
	}
}
