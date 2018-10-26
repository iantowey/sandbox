package com.lucidlogic.horsetracker.view.impl;

import android.databinding.DataBindingUtil;
import android.databinding.ViewDataBinding;
import android.os.Bundle;
import android.support.annotation.Nullable;
import android.support.v4.app.DialogFragment;
import android.view.LayoutInflater;
import android.view.View;
import android.view.ViewGroup;

import com.lucidlogic.horsetracker.BR;
import com.lucidlogic.horsetracker.R;
import com.lucidlogic.horsetracker.handler.EntityClickHandler;
import com.lucidlogic.horsetracker.handler.RaceClickHandler;
import com.lucidlogic.horsetracker.model.Entity;

/**
 * Created by ian on 01/08/16.
 */
public class EntitySelectedDialogFragment extends DialogFragment {

    public static final String LAYOUT_ID = "LAYOUT_ID";
    public static final String ENTITY = "ENTITY";

    private int layoutId;
    private Entity entity;

    public EntitySelectedDialogFragment(){};

    public static EntitySelectedDialogFragment newInstance(int layOutId, Entity entity) {
        Bundle bundle = new Bundle();
        bundle.putInt(LAYOUT_ID, layOutId);
        bundle.putParcelable(ENTITY, entity);
        EntitySelectedDialogFragment entitySelectedDialogFragment = new EntitySelectedDialogFragment();
        entitySelectedDialogFragment.setArguments(bundle);
        return entitySelectedDialogFragment;
    }

    @Override
    public void onCreate(@Nullable Bundle savedInstanceState) {
        super.onCreate(savedInstanceState);
        this.layoutId = this.getArguments().getInt(LAYOUT_ID);
        this.entity = this.getArguments().getParcelable(ENTITY);
    }

    @Override
    public View onCreateView(LayoutInflater inflater, ViewGroup container, Bundle savedInstanceState) {

        ViewDataBinding viewDataBinding = DataBindingUtil.inflate(inflater, this.layoutId, container, false);
        viewDataBinding.setVariable(BR.entity,  entity);
        viewDataBinding.setVariable(BR.handler,  new EntityClickHandler());

        return viewDataBinding.getRoot();

    }

}
